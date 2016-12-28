#!/usr/bin/env perl
use strict;
use warnings;
###########################################################################
#
#  prop_report:  A program to report information about the
#  definition of a property, including all the evidence supporting
#  the various steps.  The program also will report on the state of the
#  property for a particular database as well as the feat_names supporting
#  the various steps and those steps which are not found where appropriate.
#
#  Data is pulled from the step_feat_link, step_ev_link, prop_step
#  and prop_def tables in the common database, as well as the db_data,
#  asm_feature and ident tables of the omnium database (to get com_names).
#
#  The user must supply as an argument either a unique keyword for the
#  property name or the prop_def_id; the optional inclusion of a db name
#  as a second parameter triggers the reporting of the information for
#  that db.
#
#  The keyword "niche" is reserved, and results in output limited to the
#  biological niche ontology.
#
#  Jeremy Selengut, TIGR, 2003
#
###########################################################################

use DBI;
use Getopt::Std;
use Carp qw/confess/;
BEGIN{$SIG{__DIE__} = sub {confess @_};}

use FindBin qw/$RealBin/;

#my $server = "$RealBin/../data/common";
my $server = "common";
#my $dbtype = "SQLite";
my $dbtype = "mysql";

use constant NULL => '';

my %arg;
&getopt('f:p:e:i:a:k:n', \%arg);
my $dbfile = $arg{f};
my $prop = $arg{p} ? $arg{p} : 'all';
my $ev = $arg{e};
my $prop_def_id = $arg{i};
my $prop_acc = $arg{a};
my $keyword = $arg{k};
my $dbreport = $arg{f} ? 1 : 0;
my $niche_report = $arg{n} ? 1 : 0;
our %EVLOOKUP;

my $DEBUG = "OFF";
# my $DEBUG = "ON";

my $host = $ENV{DBSERVER} ? $ENV{DBSERVER} : 'localhost';

## login ##

my $dbh = ( DBI->connect( "dbi:$dbtype:host=$host;dbname=$server", "access", "access" ) );
if ( !defined $dbh ) {
    die "Cannot connect to server: $DBI::errstr\n";
}

## process user inputs ##
if (! $dbfile && ! $ev && ! ($prop_def_id || $prop_acc || $keyword)) {
    die "\nYou must include a dbfile (-f) and/or prop_def_id (-i), property accession (-a), keyword (-k) or HMM accession (-e).\n\n" 
}
elsif ( ! $dbfile && $ev =~ /^(TIGR|PF)\d{5}/) {
    my $query = qq{
        select distinct d.prop_def_id, d.prop_acc, d.property, ps.step_num, ps.step_name, ps.in_rule, se.get_GO 
              from prop_def d, prop_step ps, step_ev_link se
              where d.prop_def_id = ps.prop_def_id
                and ps.prop_step_id = se.prop_step_id
                and se.query = ?
    };
    my $sth = $dbh->prepare($query);
    $sth->execute($ev);
    my $results = $sth->fetchall_arrayref();
    $sth->finish;
    undef $sth;
    # print the list of properties and exit
    if ( @$results ) {
        print "\n";
        print "HMM $ev is used as evidence in the following Genome Properties:\n";
        print "\nid\taccession\tproperty\t\t\t\t\tstepnum\tstep_name\t\tin_rule\tget_GO\n";
        print "------- --------------- ----------------------------------------------- ------- ----------------------- ------- -------\n";
        foreach my $result (@$results) {
            printf "%s\t%s\t%-47.47s\t%s\t%-20.20s\t%s\t%s\n", @$result;
        }
        print "\n";
        $dbh->disconnect;
        exit(0);
    }
} elsif  ($dbfile) {
    open my $in, $dbfile || die "Can't open $dbfile for read : $!\n";
    while (my $line = <$in>) {
	chomp $line;
	my @data = split/\t/,$line;
	my $descr = pop @data;
	my @kv = split/\;/, $descr;
	foreach my $pair (@kv) {
	    my ($key, $val) = split/\=/, $pair,2;
	    if ($key eq "ID") {
		my ($ev, $acc) = split(/\_/, $val, 2);
		$EVLOOKUP{$ev}->{$data[0]} = 1;
	    }
	}
    }
}
	

my $property;
my $prop_type;
my @results;
# grab properties from db based on either accession or id or property keyword or set niche report
if ( $prop_acc =~ /^genprop/i ) { # Property Accession
    my $query = 'select prop_def_id, prop_acc, property, prop_type
                  from prop_def
                  where lower(prop_acc) = "' . $prop_acc . '"';
    @results = &do_sql( $dbh, $query );
    my $numlines = @results;
    if ( $numlines == 0 ) {
	die "\nThat property accession does not match a property.\n\n";
    }
}
elsif ( $prop_def_id =~ /^\d+/ ) { # prop_def_id
    my $query = 'select prop_def_id, prop_acc, property, prop_type 
                  from prop_def
                  where prop_def_id = ' . $prop_def_id;
    @results = &do_sql( $dbh, $query );
    my $numlines = @results;
    if ( $numlines == 0 ) {
	die "\nThat prop_def_id does not match a property.\n\n";
    }
}
elsif ($keyword) {
    my $query = 'select prop_def_id, prop_acc, property, prop_type, ispublic
                  from prop_def
                  where lower(property) like "%' . $keyword . '%"';
    @results = &do_sql( $dbh, $query );
    my $numlines = @results;
    if ( $numlines == 0 ) {
	die "\nYour keyword did not match the name of any property.\n\n";
    }
    if ( $numlines >= 2 ) {
	print "\nYour keyword was not unique and found these property names:\n\n";
	for ( my $i = 0; $i < $numlines; $i++ ) {
	    my @fields = split /\t/, $results[$i];
	    print "$fields[0]\t$fields[3]\t$fields[4]\t$fields[1]\t$fields[2]\n";
	}
	print "\n\n";
	die "Try using one of the associated property ids or property accessions\n\n";
    }
}
my @fields      = split /\t/, $results[0];
$prop_def_id = $fields[0];
$prop_acc    = $fields[1];
$property    = $fields[2];
$prop_type   = $fields[3];


## Output general property information when only a property is specified ##
if ( $dbreport == 0 || $dbreport == 1 ) {
  PROP_HEADER: {
      print "\n";
      print "-----------------------------------------------------------------------------------\n";
      print "id      accession       type    property name\n";
      print "------- --------------- ------- ---------------------------------------------------\n";
      print "$prop_def_id\t$prop_acc\t$prop_type\t$property";
      print "\n\n";
    }
  GO_LINKS: {
      my $query = qq{
	  select distinct pg.go_id, pg.get_GO
	      from prop_go_link pg
	      where pg.prop_def_id = ?
	      order by pg.get_GO
	  };
      my $sth = $dbh->prepare($query);
      $sth->execute($prop_def_id);
      my $results = $sth->fetchall_arrayref();
      if ( @$results ) {
	  print "-----------------------\n";
	  print "GO_id            get_GO\n";
	  print "---------------  ------\n";
	  foreach my $result (@$results) {
	      printf "%-16s %6s\n", @$result;
	      print "\n";
	  }
	  print "\n\n";
      }
    }
    
    ## Output general property information when no db is specified ##
    if ( 1 ) { # no db
      PROP_DEFINITION: {
	  print "-----------------------------------------------------------------------------------\n";
	  print "property definition:\n";
	  print "-----------------------------------------------------------------------------------\n";
	  my $query = 'select description from prop_def 
                      where prop_def_id = ' . $prop_def_id;
	  my @results = &do_sql( $dbh, $query );
	  &neat_print( $results[0], 80 );
	}
        
      PROP_REFS_LIT: {
	  my $query = 'select ref_num, authors, citation, title, accession, comment
                      from prop_ref
                      where prop_def_id = ' . $prop_def_id . '
                        and ref_type = "LIT"';
	  my @results = &do_sql( $dbh, $query );
	  my $numlines = @results;
	  
	  for ( my $i = 0; $i <= $numlines - 1; $i++ ) {
	      my @fields = split /\t/, $results[$i];
	      print "[$fields[0]]\t$fields[1]\n\t$fields[2]\n";
	      &neat_print_indent( $fields[3], 75 );
	      print "\t$fields[4]\n";
	      if ( ( ! $fields[5] ) or ( $fields[5] eq NULL ) ) {
		  print "\n";
	      }
	      else {
		  print "\t$fields[5]\n\n";
	      }
	  }
	}
        
      PROP_REFS_WEB: {
	  my $query = 'select url
                      from prop_ref
                      where prop_def_id = ' . $prop_def_id . '
                        and ref_type = "WEB"';
	  my @results = &do_sql( $dbh, $query );
	  my $numlines = @results;
	  unless ( $numlines == 0 ) {
	      print "-----------------------------------------------------------------------------------\n";
	      print "web links:\n";
	      print "-----------------------------------------------------------------------------------\n";
	      for ( my $i = 0; $i <= $numlines - 1; $i++ ) {
		  print "$results[$i]\n";
	      }
	      print "\n";
	  }
	}
        
      PARENT_PROP_LINKS: {
	  print "-----------------------------------------------------------------------------------\n";
	  print "parent properties                                               id      accession\n";
	  print "--------------------------------------------------------------- ------- -----------\n";
	  
	  my $query = qq{
	      select d.property, d.prop_def_id, d.prop_acc
		  from prop_def d, prop_link l
		  where l.child_id = ?
		  and d.prop_def_id = l.parent_id
		  and l.link_type = 'DAG'
	      };
	  my $sth = $dbh->prepare($query);
	  $sth->execute($prop_def_id);
	  my $results = $sth->fetchall_arrayref;
	  if ( $prop_type eq "ROOT" ) {
	      print "This property is the root of the genome properties tree.\n\n";
	  }
	  else {
	      foreach my $result (@$results) {
		  printf "%-60.60s\t%s\t%s\n", @$result;
	      }
	      print "\n";
	  }
	}

      CHILD_PROP_LINKS: {
	  my $query = qq{
	      select d.property, d.prop_def_id, d.prop_acc
		  from prop_def d, prop_link l
		  where l.parent_id = ?
		  and d.prop_def_id = l.child_id
		  and l.link_type = 'DAG'
	      };
	  my $sth = $dbh->prepare($query);
	  $sth->execute($prop_def_id);
	  my $results = $sth->fetchall_arrayref();
	  if ( @$results ) {
	      print "-----------------------------------------------------------------------------------\n";
	      print "child properties                                                id      accession\n";
	      print "--------------------------------------------------------------- ------- -----------\n";
	      foreach my $result (@$results) {
		  printf "%-60.60s\t%s\t%s\n", @$result;
	      }
	      print "\n";
	  }
	}
        
      IMPLIED_PROP_LINKS: {
	  my $query = qq{
	      select d.property, 'NO', d.prop_def_id, d.prop_acc
		  from prop_def d, prop_link l
		  where l.parent_id = ?
		  and d.prop_def_id = l.child_id
		  and l.link_type = 'IMPLIED_BY'
		  union
		  select d.property, 'YES', d.prop_def_id, d.prop_acc
		  from prop_def d, prop_link l
		  where l.child_id = ?
		  and d.prop_def_id = l.parent_id
		  and l.link_type = 'IMPLIED_BY'
	      };
	  my $sth = $dbh->prepare($query);
	  $sth->execute($prop_def_id, $prop_def_id);
	  my $results = $sth->fetchall_arrayref;
	  if (@$results) {
	      print "-----------------------------------------------------------------------------------\n";
	      print "properties to which state is applied by inference       state   id      accession\n";
	      print " (and from which the opposite state is inferred)\n";
	      print "------------------------------------------------------- ------- ------- -----------\n";
	      foreach my $result (@$results) {
		  printf "%-50.50s\t%s\t%s\t%s\n", @$result;
	      }
	      print "\n";
	  }
	}

      UPSTREAM_PROP_LINKS: {
	  my $query = qq{
	      select d.property, d.prop_def_id, d.prop_acc
		  from prop_def d, prop_link l
		  where l.child_id = ?
		  and d.prop_def_id = l.parent_id
		  and l.link_type = "BRANCH_POINT"
	      };
	  my $sth = $dbh->prepare($query);
	  $sth->execute($prop_def_id);
	  my $results = $sth->fetchall_arrayref();
	  if (@$results) {
	      print "-----------------------------------------------------------------------------------\n";
	      print "upstream properties                                             id      accession\n";
	      print "--------------------------------------------------------------- ------- -----------\n";
	      foreach my $result (@$results) {
		  printf "%-60.60s\t%s\t%s\n", @$result;
	      }
	      print "\n";
	  }
	}
        
      DOWNSTREAM_PROP_LINKS: {
	  my $query = qq{
	      select d.property, d.prop_def_id, d.prop_acc
		  from prop_def d, prop_link l
		  where l.parent_id = ?
		  and d.prop_def_id = l.child_id
		  and l.link_type = 'BRANCH_POINT'
	      };
	  my $sth = $dbh->prepare($query);
	  $sth->execute($prop_def_id);
	  my $results = $sth->fetchall_arrayref;
	  
	  if (@$results) {
	      print "-----------------------------------------------------------------------------------\n";
	      print "downstream properties                                           id      accession\n";
	      print "--------------------------------------------------------------- ------- -----------\n";
	      foreach my $result (@$results) {
		  printf "%-60.60s\t%s\t%s\n", @$result;
	      }
	      print "\n";
	  }
	}
        
      PROP_STEPS: {
	  my $step_q = 'select 1 from prop_step where prop_def_id = ' . $prop_def_id;
	  my $has_steps = &do_sql( $dbh, $step_q );
	  if ($has_steps) {
	      print "-----------------------------------------------------------------------------------\n";
	      print "stepnum in_rule step_name\n";
	      print "------- ------- -------------------------------------------------------------------\n";
	      
	      my $step_query = qq{
		  select step_num, in_rule, step_name
		      from prop_step
		      where prop_def_id = ?
		      order by step_num
		  };
	      my $sth = $dbh->prepare($step_query);
	      $sth->execute($prop_def_id);
	      my $step_results = $sth->fetchall_arrayref({});
	      foreach my $result (@$step_results) {
		  printf "%-7.7s\t%1.1s\t%-.60s\n", @$result{qw/step_num in_rule step_name/};
	      }
	      print "\n";
	      
	      print "-----------------------------------------------------------------------------------\n";
	      print "step_num get_GO method          query\n";
	      print "-------- ------ --------------- -----------------------\n";
	      
	      foreach my $step (@$step_results) {
		  my $step_num = $step->{step_num};
		  my $step_name = $step->{step_name};
		  print "$step_num\t($step_name)\n";
		  
		  my $query = qq{
		      select se.get_GO, se.method, se.query
			  from step_ev_link se, prop_step ps
			  where ps.step_num = ?
			  and ps.prop_def_id = ?
			  and ps.prop_step_id = se.prop_step_id
			  and se.query not like '%CLUSTER%'
			  and se.query not like '%BLAST%'
			  and se.query not like '%FUSION%'
			  and se.query not like '%ALIGN%'
		      };
		  my $sth = $dbh->prepare($query);
		  $sth->execute($step_num, $prop_def_id);
		  my $results = $sth->fetchall_arrayref();
		  
		  foreach my $result ( @$results) {
		      printf "         %-6.6s\t%-8.8s\t%-23.23s", @$result;
		      if ($dbfile) { &print_hits(\%EVLOOKUP, @$result) }
		      print "\n";
		  }
		  print "\n";
	      }
	  }
	  else {
	      print "\nThis property is of type " . $prop_type . " and has no steps to report.\n\n";
	  }
	}
    }

    ## Output property information specific to the input db and property ##

    else {
        # TODO Get results for an individual property/db working
#       CHECK_FOR_DB: {
#           my $query   = 'select 1 from property where db = "' . $key2 . '"';
#           my @results = &do_sql( $dbh, $query );
#           my $dbfound = @results;
#           if ( $dbfound == 0 ) {
#               die "That database name does not exist in the genome properties tables,
#   check the spelling and try again.\n\n";
#           }
#       }
#       
#       CHECK_FOR_PROP_DEF_ID: {
#           my $query      = 'select 1 from property where prop_def_id = ' . $prop_def_id;
#           my @results    = &do_sql( $dbh, $query );
#           my $prop_found = @results;
#           if ( $prop_found == 0 ) {
#               die "That property does not have assigned states.\n\n";
#           }
#       }
#       
#       SHOW_PROPERTY_STATE: {
#           my $prop_state_query = qq{
#               select p.state, p.value, p.assignby 
#                     from property p, prop_def d
#                     where p.property = d.property
#                       and d.prop_def_id = ?
#                       and p.db = ?
#           };
#           my $prop_sth = $dbh->prepare($prop_state_query);
#           $prop_sth->execute($prop_def_id, $key2);
#           my @prop_state = $prop_sth->fetchrow_array();
#           $prop_sth->finish;
#           undef $prop_sth;
#           if (! @prop_state) {
#               print "This property has not been evaluated for this genome.\n\n";
#           }
#           else {
#               print "-----------------------------------------------------------------------------------\n";
#               print "state                                   value           assignby\n";
#               print "--------------------------------------- --------------- ---------------------------\n";
#               printf "%-35.35s\t%s\t\t%s\n", @prop_state;
#   
#               SHOW_COMMENTS: {
#                   my $query = 'select p.auto_comment, p.curator_comment
#                         from property p, prop_def d
#                         where p.property = d.property
#                           and d.prop_def_id = ' . $prop_def_id . '
#                           and p.db = "' . $key2 . '"';
#   
#                   my @results = &do_sql( $dbh, $query );
#                   my ( $auto_comment, $curator_comment ) = split /\t/, $results[0];
#                   if ( $auto_comment ne "" ) {
#                       print "-----------------------------------------------------------------------------------\n";
#                       print "RULE-GENERATED COMMENT:\n";
#                       print "-----------------------------------------------------------------------------------\n";
#                       &neat_print( "$auto_comment\n", 80 );
#                   }
#                   chomp $curator_comment;
#                   if ( $curator_comment ne "" ) {
#                       print "-----------------------------------------------------------------------------------\n";
#                       print "CURATOR COMMENT:\n";
#                       print "-----------------------------------------------------------------------------------\n";
#                       &neat_print( $curator_comment, 80 );
#                   }
#               }
#               
#               SHOW_STEPS: {
#                   if ( ( $prop_type eq "PATHWAY" ) or ( $prop_type eq "SYSTEM" ) or ( $prop_type eq "GUILD" ) or ( $prop_type eq "TEST" ) or ( $prop_type eq "PATHEMA" ) or ( $prop_type eq "METAPATH" ) or ( $prop_type eq "COUNTER" ) or ( $prop_type eq "OPERON" ) ) {
#       
#                       my $step_num_query = qq{
#                           select distinct step_num, step_name
#                                 from prop_step 
#                                 where prop_def_id = ' . $prop_def_id . ' 
#                                 order by step_num
#                       };
#                       my $step_num_sth = $dbh->prepare($step_num_query);
#                       my $step_nums = $step_num_sth->execute($prop_def_id);
#                       my $step_nums = $step_num_sth->fetchall_arrayref({});
#                       print "---------------------------------------------------------------------------------------------------------------------------------\n";
#                       print "stepnum feat_name       locus         query     method        assignby      com_name\n";
#                       print "------- --------------- ------------- --------- ------------- ------------- -----------------------------------------------------\n";
#                       my $missing = 0;
#                       my @missing;
#                       foreach my $step (@$step_nums) {
#                           my $step_num = $step->{step_num};
#       
#                           my $standard_query = 'select distinct convert(char(12), sf.feat_name), 
#                                                convert(char(12), dbi.locus) + "  " +
#                                                convert(char(15), se.query),
#                                                convert(char(10), se.method),
#                                                convert(char(10), sf.assignby), 
#                                                convert(char(42), dbi.com_name)
#                                         from step_feat_link sf, step_ev_link se, prop_step ps, 
#                                              ' . $key2 . '..asm_feature dba, ' . $key2 . '..stan dbs, ' . $key2 . '..ident dbi
#                                         where sf.step_ev_id = se.step_ev_id
#                                           and se.prop_step_id = ps.prop_step_id
#                                           and ps.prop_def_id = ' . $prop_def_id . '
#                                           and ps.step_num = "' . $step_nums[$i] . '"
#                                           and sf.db = "' . $key2 . '" 
#                                           and sf.feat_name = dba.feat_name
#                                           and dba.feat_name = dbi.feat_name
#                                           and dba.asmbl_id = dbs.asmbl_id
#                                           and dbs.iscurrent = 1
#                                         union all
#                                         select distinct convert(char(12), sf.feat_name),  
#                                                convert(char(12), dbi.locus) + "  " +
#                                                convert(char(15), se.query),
#                                                convert(char(10), se.method),
#                                                convert(char(10), sf.assignby), 
#                                                convert(char(42), dbi.com_name)
#                                         from step_feat_link sf, step_ev_link se, prop_step ps, 
#                                              ' . $key2 . '..asm_feature dba, ' . $key2 . '..stan dbs, ' . $key2 . '..nt_ident dbi
#                                         where sf.step_ev_id = se.step_ev_id
#                                           and se.prop_step_id = ps.prop_step_id
#                                           and ps.prop_def_id = ' . $prop_def_id . '
#                                           and ps.step_num = "' . $step_nums[$i] . '"
#                                           and sf.db = "' . $key2 . '" 
#                                           and sf.feat_name = dba.feat_name
#                                           and dba.feat_name = dbi.feat_name
#                                           and dba.asmbl_id = dbs.asmbl_id
#                                           and dbs.iscurrent = 1';
#       
#                           my $trna_query = 'select se.method 
#                                     from step_ev_link se, prop_step ps
#                                     where se.method like "%RNA"
#                                       and se.prop_step_id = ps.prop_step_id
#                                       and ps.step_num = "' . $step_nums[$i] . '"';
#       
#                           my @trna_results = &do_sql( $dbh, $trna_query );
#                           my $tRNAlines = @trna_results;
#                           my @results;
#                           if ( $tRNAlines != 0 ) {
#                               my $query = 'select distinct convert(char(14), sf.feat_name), + "              " + convert(char(15), se.query), 
#                                                sf.assignby
#                                         from step_feat_link sf, step_ev_link se, prop_step ps, 
#                                              ' . $key2 . '..asm_feature dba, ' . $key2 . '..stan dbs
#                                           where sf.step_ev_id = se.step_ev_id
#                                           and se.prop_step_id = ps.prop_step_id
#                                           and ps.prop_def_id = ' . $prop_def_id . '
#                                           and ps.step_num = "' . $step_nums[$i] . '"
#                                           and sf.db = "' . $key2 . '" 
#                                           and sf.feat_name = dba.feat_name
#                                           and dba.asmbl_id = dbs.asmbl_id
#                                           and dbs.iscurrent = 1';
#       
#                               @results = &do_sql( $dbh, $query );
#                           }
#                           elsif ( $prop_type eq "METAPATH" ) {
#                               my $query = 'select convert(char(15), p.state), "--", convert(char(15), se.query), convert(char(10), p.assignby), ""
#                                         from step_ev_link se, prop_step ps, property p, prop_def d
#                                         where se.query = d.prop_acc
#                                           and d.prop_def_id = p.prop_def_id
#                                           and p.db = "' . $key2 . '"
#                                           and se.prop_step_id = ps.prop_step_id
#                                           and ps.step_num = "' . $step_nums[$i] . '"
#                                           and ps.prop_def_id = ' . $prop_def_id . '
#                                         union
#                                         ' . $standard_query;
#       
#                               #print "\n$query\n";
#                               @results = &do_sql( $dbh, $query );
#                           }
#                           else {
#       
#                               my $query = $standard_query;
#                               @results = &do_sql( $dbh, $query );
#                           }
#
#                           my $newnumlines = @results;
#                           if ( $newnumlines == 0 ) {
#                               $missing[$missing] = $step;
#                               $missing++;
#                           }
#                           else {
#                               my $foundone = "FALSE";
#                               foreach my $line (@results) {
#                                   unless ( ( $line =~ /NOCLUST/ ) && ( $line =~ /manual/ ) ) {
#                                       $foundone = "TRUE";
#                                   }
#                               }
#                               if ( $foundone eq "FALSE" ) {
#                                   $missing[$missing] = $step;
#                                   $missing++;
#                               }
#                               else {
#                                   print $step;
#                                   for ( my $j = 0; $j <= $newnumlines - 1; $j++ ) {
#                                       print "\n";
#                                       print "        " . $results[$j];
#                                   }
#                                   print "\n\n";
#                               }
#                           }
#                       }
#                       print "\n";
#                       if ( $missing != 0 ) {
#                           print "-----------------------------------------------------------------------------------\n";
#                           print "These steps are missing:\n\n";
#                           print "in_rule stepnum step name\n";
#                           print "------- ------- -------------------------------------------------------------------\n";
#                           for ( my $i = 0; $i <= $missing - 1; $i++ ) {
#                               my $step_num = $missing[$i];
#                               $step_num =~ s/\t.+//;
#                               my $query   = 'select in_rule from prop_step where prop_def_id = ' . $prop_def_id . ' and step_num like "' . $step_num . '%"';
#                               my @results = &do_sql( $dbh, $query );
#                               my $in_rule = $results[0];
#                               $missing[$i] =~ s/\t//;
#                               $missing[$i] =~ s/\(|\)//g;
#                               print $in_rule. "\t" . $missing[$i] . "\n";
#                               my @parts           = split /\s/, $missing[$i];
#                               my $missed_step_num = $parts[0];
#                               $query           = 'select se.query
#                                         from step_feat_link sf, step_ev_link se, prop_step ps
#                                         where ps.step_num = "' . $missed_step_num . '"
#                                           and ps.prop_step_id = se.prop_step_id
#                                           and se.step_ev_id = sf.step_ev_id
#                                           and sf.db = "' . $key2 . '"
#                                           and sf.feat_name = "MISSING"
#                                           and ps.prop_def_id = ' . $prop_def_id;
#                               ##print "\n$query\n";
#                               @results = &do_sql( $dbh, $query );
#                               foreach my $result (@results) {
#                                   print "\t\t$result checked by allsearch and not found\n";
#                               }
#                           }
#                           print "\n";
#                       }
#                       my $cluster_query = 'select max(cluster_num) 
#                                         from cluster_link 
#                                         where db = "' . $key2 . '"
#                                           and prop_def_id = ' . $prop_def_id;
#                       my @results = &do_sql( $dbh, $cluster_query );
#                       my $number_of_clusters;
#                       if   ( $results[0] eq NULL ) { 
#                           $number_of_clusters = 0 
#                       }
#                       else { 
#                           $number_of_clusters = $results[0] 
#                       }
#                       print "These features are organized into $number_of_clusters multi-gene co-localized clusters.\n\n";
#                   }
#                   else { 
#                       print "\nThis property is of type " . $prop_type . " and has no steps to report.\n\n"; 
#                   }
#               }
#           }
#       }
    }
}

## Output a report of properties for a specified genome ordered by the DAG ##
else {
    my $root_id;
    my $restrict = 1;
    #print "/n$key2/n";
    if ( $niche_report ) {
        $root_id        = 4005;    ## biological niche
        $restrict = 0;
    }
    else {
        if   ( uc($prop_def_id) eq "ALL" ) { 
            $restrict = 0 
	}
        my $sth = $dbh->prepare("SELECT prop_def_id from prop_def where prop_type = 'ROOT'");
        $sth->execute;
        ($root_id) = $sth->fetchrow_array();
        $sth->finish;
    }
    
    print_node_tree($dbh, $root_id, $dbfile, {restrict=>$restrict});

}

## disconnect ##
$dbh->disconnect;
exit(0);

###subroutines ###

sub print_node_tree {
    my $dbh = shift;
    my $root_id = shift;
    my $db = shift;
    my $options = shift;
    my $return_root = $options->{return_root};
    my $level = $options->{level} || 0;
    my $restrict = $options->{restrict} || 1;
    
    my @results;
    my @children;
    if ($return_root) {
	@children =  $root_id;
    }
    else {
	@children = _get_child_nodes($dbh, $root_id);
    }
    foreach my $node (@children) {
	my $child_info = _get_node_info($dbh, $node, $db);
	if (!$restrict or $child_info->{ispublic}) {
	    $child_info->{level} = $level; 
	    if ($level == 0) {
		$child_info->{name} = uc($child_info->{name});
	    }
	    my $indent = '    ' x $child_info->{level};
	    my $prop_def_id = $child_info->{prop_def_id};
	    my $prop_name = $indent.($child_info->{name}||'');
	    my $prop_type = $child_info->{type}||'';
	    my $prop_state = $child_info->{state}||'';
	    my $prop_value = $child_info->{value}||'';
	    if ($prop_value) { $prop_value = sprintf('%1.4f', $prop_value); };
	    my $assignby = $child_info->{assignby} || '';
	    printf "%7s\t%-56.56s\t%6.6s\t%10.10s\t%6.6s\t%s\n", $prop_def_id, $prop_name, $prop_type, $prop_state, $prop_value, $assignby;
	    push @results, $child_info; 
	    push @results, print_node_tree($dbh, $node, $db, {level=>($level+1), $restrict=>$restrict});
	}
    }
    return @results;
}

sub _get_child_nodes {
    my $dbh = shift;
    my $prop_def_id = shift;
    my $q = qq{
        select l.child_id
	    from prop_link l
	    where l.link_type = 'DAG' 
	    AND l.parent_id = ? 
	};
    my $sth = $dbh->prepare($q);
    $sth->execute($prop_def_id);
    my @children = map {$_->[0]} @{$sth->fetchall_arrayref()};
    return @children;
}

sub _get_node_info {
    my $dbh = shift;
    my $prop_def_id = shift;
    my $db = shift;
    my $q = qq{
        select d.prop_def_id, d.property as name, d.prop_type as type, d.ispublic, 
	p.state, p.value, p.assignby
	    from prop_def d
	    left join ( 
			select property, state, value, assignby
			from property
			where db = ?
			) as p
			on d.property = p.property
			where d.prop_def_id = ?
		    };
    my $sth = $dbh->prepare($q);
    $sth->execute($db, $prop_def_id);
    my $info = $sth->fetchrow_hashref();
    $sth->finish;
    return $info; 
}

## do a search
sub do_sql {
    my ( $dbh, $query, @params ) = @_;
    my ( @results, @row );

    if ( $DEBUG eq "ON" ) {
        print "\n\n$query\n\n";
    }
    my $sth = $dbh->prepare($query);
    if ( !defined $sth ) {
        die "Cannot prepare statement: $DBI::errstr\n";
    }
    $sth->execute(@params) || die $sth->errstr;

    #Fetch the rows back from the SELECT statement;
    while ( @row = $sth->fetchrow() ) {
        push( @results, join( "\t", map {defined($_)?$_:NULL} @row ) );
    }

    # Release the statement handle resources;
    $sth->finish;

    if ( $DEBUG eq "ON" ) {
        print "@results\n\n";
        my $reply = <STDIN>;
    }

    return @results;
}

sub wrap_line {
    my $line = shift;
    my $max_length = shift;
    my $paragraph = $line;
    $paragraph =~ s/([^\n]{50,70})(?:\b\s*|\n)/$1\n/gi;
    while ( $paragraph =~ m/([^\n]{71})/ ) {
        $paragraph =~ s/([^\n]{70})([^\n])/$1\n$2/g;
    }
    return $paragraph;
}

sub neat_print {
    my ( $string, $length ) = @_;
    my @lines = split /\n/, $string;
    foreach my $line (@lines) {
        if (length($line) > $length) {
            $line = wrap_line($line, $length);
        }
        print "$line\n";
    }
}

sub neat_print_indent {
    my ( $string, $length ) = @_;
    my @lines = split /\n/, $string;
    foreach my $line (@lines) {
        $line = wrap_line($line, $length);
        $line =~ s/\n/\n\t/g;
        print "\t$line\n";
    }
}

sub print_hits {
    my $EVLOOKUP = shift;
    my $ev = pop @_;
    if (defined $EVLOOKUP->{$ev}) {
	print "\t", join(";", keys %{$EVLOOKUP->{$ev}});
    }
}

=pod

    =head2 print_node_tree
    given a root prop_def_id, follows the dag and returns the tree
    
    =head3 parameters

  1: dbh = database handle
  2: root_id = prop_def_id of root node
  3: db = database name to return results for
  4: options = hashref {
    return_root => bool, true if root should also be returned, defaults to false
    level => number, defaults to 0,
    restrict => bool, true if restricted to public entries, defaults to true

    =head3 returns (in expected order of display)

    (
{
    level   => <dist_from_root>, 
    prop_def_id => <prop_def_id>, 
    name    => <name of property>,
    type    => <property type>,
    state   => <property state (undef if unevaluable)>
	value   => <property value (undef if unevaluable)>
	assignby => <property assigned by (undef if unevaluable)>
},
    ...
    )

    =cut
    
