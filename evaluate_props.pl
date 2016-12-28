#!/usr/bin/env perl
use warnings;
use strict;
use Carp;


use version; our $VERSION = qv( qw$Revision 0.2.1$[1] );
use Log::Log4perl qw(:levels);

use Getopt::Long qw(:config no_ignore_case); # Option parser. Exports GetOptions()
use Pod::Usage; # for --help parameter. Exports pod2usage(), which exit()s
# use IO::Prompt; # for command input. Exports prompt()

## Uncomment if script uses lib directory below script directory
use FindBin;
use lib ($FindBin::RealBin."/../blib/lib",$FindBin::RealBin."/../lib");

# Local Modules
use GenProp;
use GenProp::DB;
use GenProp::Exception;
use GenProp::Genome::Factory;
use GenProp::Genome::GFF;
use GenProp::Property::Factory;
use GenProp::Property::Evaluated::Factory;

sub read_file {
    # What file are we reading here?
    my $filename = shift;
    open my $fh, "<", $filename or die "Can't open $filename to read. $!\n";
    my @lines = grep {$_} map { s/^\s+//; s/#.*$//; s/\s+$//; $_ } <$fh>;
    close $fh or die "Can't close $filename. $!\n";
    return @lines;
}

sub parse_options {

    my $opts = {};
    $opts->{update_db} = 1; # default to record
    my @prop_or_rule_ids;
    my @genomes;
    GetOptions(
        'genome_db|d=s'         => sub {push @genomes, $_[1]},
        'file_genome_db|fd=s'   => sub {push @genomes, read_file($_[1])},

        'rule_id|i=i'           => sub {push @prop_or_rule_ids, {id => $_[1], type => 'rule'}; },
        'rule_type|T=s'         => sub {push @prop_or_rule_ids, {id => $_[1], type => 'rule_type'};},
        'prop_id|p=i'           => sub {push @prop_or_rule_ids, {id => $_[1], type => 'prop'};},
        'accession|a=s'         => sub {push @prop_or_rule_ids, {id => $_[1], type => 'accession'};},
        'step_id|s=i'           => sub {push @prop_or_rule_ids, {id => $_[1], type => 'step'};},
        'file_prop_id|fp=s'        => sub {push @prop_or_rule_ids, map { {id=>$_, type => 'prop'} } read_file($_[1])},
        'file_accession|fa=s'   => sub {push @prop_or_rule_ids, map { {id=>$_, type => 'accession'} } read_file($_[1])},

        'config|C=s'            => sub {$opts->{config_file} = $_[1] },
        'username|U=s'          => sub {warn "--username is ignored. Use --config to set config file\n";},
        'password|P=s'          => sub {warn "--password is ignored. Use --config to set config file\n";},
        'server|S=s'            => sub {warn "--server is ignored. Use --config to set config file\n";},

        'verbose|debug|D+'      => sub {$opts->{verbose} = $_[1] },
        'quiet|q'               => sub { $opts->{'quiet'} = $_[1] },

        'execute|E!'            => sub {$opts->{update_db} = $_[1] },
        'dry_run'               => sub { $opts->{update_db} = 0 },

        'follow'                => sub { $opts->{follow} = 1 },

        'help|?'                => sub {die pod2usage( -exitstatus => 0)},
        'version'               => sub {print "$0 version $VERSION\n"; exit 0;},
        'man'                   => sub {die pod2usage( -exitstatus => 0, -verbose => 1)},
    ) or pod2usage(2);
    $opts->{props} = \@prop_or_rule_ids;
    $opts->{genomes} = \@genomes;
    return $opts;
}

sub get_all_evaluable_prop_ids {
    ## INFO this is a pretty nasty hack until I get property navigation working property
    my $db = GenProp::DB->new('GenPropDef');
    my $dbh = $db->get_dbh;
    my $sql = qq(
        SELECT prop_def_id
        FROM prop_def
        WHERE eval_method IN ('RULES', 'SIMPLE','GUILD_SUFF')
        AND ispublic != -1
        ORDER BY prop_def_id
    );
    my $res = $dbh->selectall_arrayref($sql);
    my @prop_ids;
    if (@$res) {
        @prop_ids = map {$_->[0]} @$res;
    }
    return @prop_ids;
}

sub get_prop_id_from_rule_id {
    my $rule_id = shift;
    my $db = GenProp::DB->new('GenPropDef');
    my $dbh = $db->get_dbh;
    my $sql = qq(
        SELECT prop_def_id
        FROM rules
        WHERE rule_type like 'SGC_%'
        AND id = ?
    );
    my $sth = $dbh->prepare($sql);
    $sth->execute($rule_id);
    my ($prop_id) = $sth->fetchrow_array;
    $sth->finish;
    return $prop_id;
}

## Main script implementation here


# parse command-line options
my $options_r = parse_options();
    # this should result in
    # a genome id
    # 1+ property/rule ids to evaluate
    # flags: verbose, quiet, update_db
    # optionally, a config file
    # could also exit if documentation requested

if ( ! $options_r->{genomes}) {
    die "No genome passed. Cannot Evaluate\n";
}


# initialize app (load config info, etc)
my $app = GenProp->init($options_r->{config_file});
my $logger = GenProp->get_logger;
if ($options_r->{verbose}) {
    $logger->level($DEBUG);
}

$logger->debug('Retrieving Properties');
my @prop_ids;
# expand all prop_ids request
if (@{$options_r->{props}} == 1 and $options_r->{props}[0]{id} eq 'ALL') {
    shift @{$options_r->{props}}; # dump the prop_id 'ALL'
    @prop_ids = map{ {id => $_, type => 'prop'} } get_all_evaluable_prop_ids();
}
# convert rule ids and accessions to property ids
foreach my $id_r ( @{ $options_r->{props} } ) {
    if ($id_r->{type} eq 'prop' or $id_r->{type} eq 'accession' or $id_r->{type} eq 'step' ) {
        push @prop_ids, $id_r;
    }
    elsif ($id_r->{type} eq 'rule') {
        my $prop_id = get_prop_id_from_rule_id($id_r->{id});
        push @prop_ids, {id => $prop_id, type => 'prop'};
    }
    elsif ($id_r->{type} eq 'rule_type') {
        my $db = GenProp::DB->new('GenPropDef');
        my $dbh = $db->get_dbh;
        my $sth = $dbh->prepare(qq{
            SELECT prop_def_id
            FROM rules
            WHERE rule_type = ?
            ORDER BY prop_def_id
        });
        $sth->execute($id_r->{id});
        while (my ($prop_id) = $sth->fetchrow_array) {
            push @prop_ids, {id => $prop_id, type => 'prop'};
        }
    }
    else {
        croak "Don't know what type (",$id_r->{type},") is.\nThis is a code error.";
    }
}

# expand genomes list
if ( @{$options_r->{genomes}} == 1 and $options_r->{genomes}[0] eq 'ALL') {
    @{$options_r->{genomes}} = GenProp::Genome::Factory->get_genome_db_names({private=>1});
}


if ( ! @{$options_r->{genomes}} ) {
    die "No genomes specified. Nothing to evaluate.\n";
}
print "Database\tPropAcc\tPropName\tState\tVal\tAssigned-by\tComment\n";
DATABASE:
foreach my $db_name (@{$options_r->{genomes}}  ) {
    $logger->debug("Getting genome and model for $db_name");
    my $genome;
    my $gff_files;
    if ($db_name =~ /:/) {
        my @files = split(/:/, $db_name);
        $db_name = shift @files;
        eval {
            $genome = GenProp::Genome::GFF->load({db_name => $db_name, gff_files=>\@files});
        };
        if (my $e = $@) {
            $logger->warn("Error loading genome $db_name. $@");
            next DATABASE;
        }
    }
    else {
        eval {
            $genome = GenProp::Genome::Factory->load({ genome_id => $db_name });
        };
        my $e;
        if ($e = GenProp::Exception::Genome::InvalidName->caught) {
            $logger->warn("Genome $db_name does not exist. Skipping.");
            next DATABASE;
        }
        elsif ($e = GenProp::Exception::DB::ConnectionFailed->caught) {
            $logger->warn("Could not connect to database $db_name. Skipping\nError was: ".$e->error);
            next DATABASE;
        }
        elsif ($@) {
            die $@;
        }
    }

    $logger->info('Evaluating '.scalar(@prop_ids).' properties against '.$genome->db_name);
    foreach my $id_r (@prop_ids) {
        my $property;
        if ($id_r->{type} eq 'prop') {
            $property = GenProp::Property::Factory->load({ def_id => $id_r->{id}});
        }
        elsif ($id_r->{type} eq 'accession') {
            $property = GenProp::Property::Factory->load({ accession => $id_r->{id}});
        }
        elsif ($id_r->{type} eq 'step') {
            $property = GenProp::Property::Factory->load({ step_id => $id_r->{id} });
        }

        if ( !$property) {
            $logger->warning('Property ID $prop_id does not exist\n');
            next;
        }
        my $prop_name = $property->name;
        my $prop_id = $property->id;
        my $prop_eval_method = $property->eval_method || '';
        $logger->info("Evaluating $prop_name ($prop_id) via $prop_eval_method\n");
        my $evald_prop = GenProp::Property::Evaluated::Factory->load({property=>$property, genome=>$genome});
        eval {
            if ($options_r->{follow}) {
                $evald_prop->evaluate({eval_child_if_older_than=>$app->init_time()});
            }
            else {
                $evald_prop->evaluate();
            }

            if ($options_r->{update_db} and !($prop_eval_method eq 'RULES')) {
                $evald_prop->record;
            }
        };
        if (my $e = $@) {
            $logger->error("Could not evaluate $prop_name ($prop_id).\n$e\n----");
        }

        if ($evald_prop->state) {
            print join("\t", ($evald_prop->genome->db_name,
			      $evald_prop->accession,
			      $evald_prop->name,
			      $evald_prop->state,
			      ( defined($evald_prop->value)?$evald_prop->value:'NULL' ),
			      $evald_prop->assigned_by,
			      ( $evald_prop->auto_comment || "" ))),
	    "\n";
	    foreach my $step_ref (sort {$a->rank <=> $b->rank || $a->rank cmp $b->rank} $evald_prop->get_step_links) {
		my @genome_features;
		print "\t", $step_ref->child->rank, "\t",
		$step_ref->child->name, "\t";
		foreach my $ev($step_ref->child->get_evidence) {
		    my @features = $genome->find_feature_by_hmm($ev->{parameters}->{query});
		    if (@features) {
			print join(";", map { $_->{feat_name} } @features), ";";
#			print "($ev->{parameters}->{query});"; # right now, the feat_name above contains the hmm acc.
		    }
		}
		print "\n";
	    }
        }
    }
    $genome->uncache;
}

__END__

=head1 NAME

evaluate_props.pl - given a small genome name and a property id, rule id, or rule type,
applies that rule or rules of that type to the genome


=head1 VERSION

This document describes evaluate_props.pl version 2.0.1


=head1 SYNOPSIS

    evaluate_props.pl (-d genome_database, ...)|(-fd file_of_genome_db_names) ( -p property_id | -fp file_of_prop_ids | -i rule_id | -T rule_type | -a prop_accession | -fa file_of_prop_accessions ) [-D] [--dry_run] [--follow] [-h] [-man]
    At least one genome_db and at least one of -p, -i, or -T are required.

=head1 OPTIONS

=over 8

=item B<--config> or B<-C>

Override GenProp config file name. Defaults to searching for
GenProp.ini through ~/.genprop:../conf:./conf:. (based on location of the this program).

=item B<--genome> or B<-d>

small genome database name -- Multiple databases are OK.
Use 'ALL' to run against all genomes already loaded into Genome Properties.
To specify GFF files to use as the genome database instead of loading from
a relational database, use <<db_name>>:<<file>>[:<<file>> ...]

=item B<--file_genome_id> or B<-fd>

file of genome database IDs, one per line

=item B<--prop_id> or B<-p>

property ID to evaluate -- Multiple prop IDs are OK.

=item B<--file_prop_id> or B<-fp>

file of property IDs, one per line

=item B<--accession> or B<-a>

property accession to evaluate -- Multiple accessions are OK.
Use 'ALL' to evaluate all properties.

=item B<--file_accession> or B<-fa>

file of accessions, one per line

=item B<--step_id> or B<-s>

Prop_Step ID to check -- Multiple prop_step IDs are OK.

=item B<--rule_id> or B<-i>

rule id to run -- Multiple rule IDs are OK.

=item B<--rule_type> or B<-T>

type of rules to run all of -- Multiple types are OK.

=item B<--debug> or B<-D>

print verbose messages.

=item B<--execute> or B<-E>

update the database. Currently always on

=item B<--dry_run>

opposite of --execute. Currently never true

=item B<--follow>

evaluate all dependencies. Without --follow, will evaluate dependencies that have not been updated in the last 24 hours

=item B<--help> or B<-h>

print a brief help message and exit

=item B<--man>

print the man page for this script

=back

=head1 DESCRIPTION

evaluate_props.pl is designed to be a replacement for Dan Haft's original evaluate_props.pl.
Given a rule (or type of rule) to run and a small genome database to run
against, it will update the genome properties database with the information
for that (or those) rules.


=head1 DIAGNOSTICS

=for author to fill in:
    List every single error and warning message that the module can
    generate (even the ones that will "never happen"), with a full
    explanation of each problem, one or more likely causes, and any
    suggested remedies.

=over

=item C<< Error message here, perhaps with %s placeholders >>

[Description of error here]

=item C<< Another error message here >>

[Description of error here]

[Et cetera, et cetera]

=back


=head1 CONFIGURATION AND ENVIRONMENT

evaluate_props.pl requires GenProp.ini


=head1 DEPENDENCIES

evaluate_props.pl is dependent on the Sybase databases at TIGR, particularly
the small genome databases and the common database. The functions run from
evaluate_props.pl are dependent on the hmm tables in the egad database.

=head1 BUGS AND LIMITATIONS

=for author to fill in:
    A list of known problems with the module, together with some
    indication Whether they are likely to be fixed in an upcoming
    release. Also a list of restrictions on the features the module
    does provide: data types that cannot be handled, performance issues
    and the circumstances in which they may arise, practical
    limitations on the size of data sets, special cases that are not
    (yet) handled, etc.

=over 8

=item B<No way to do dry run>

Because most rules are running SQL statements that directly update the
genome properties tables, there is no way to do a dry run.

=item B<Default user and database crash evaluate_props.pl>

The default user and password ('access') do not have write permission
to the databases, and so the SQL statements run by the rules will fail,
and evaluate_props.pl will crash.

=back

Please report any bugs or feature requests to Alexander Richter C<< <arichter@tigr.org> >>

=head1 AUTHOR

Alexander Richter  C<< <arichter@tigr.org> >>


=head1 LICENCE AND COPYRIGHT

Copyright (c) 2006, Alexander Richter C<< <arichter@tigr.org> >> and TIGR.
All rights reserved.

This module is free software; you can redistribute it and/or
modify it under the same terms as Perl itself. See L<perlartistic>.


=head1 DISCLAIMER OF WARRANTY

