#!/usr/bin/env perl

if ( $#ARGV > -1 ){   # collect arguments (nems of files to process
 @listfile = (@ARGV) ;
 $listfile[$#listfile] =~ s/\n/ /;
}
$use_preproc=1;
foreach $file (@listfile) {
  if($file =~ /^-h|--help/) {
    print STDERR "usage : mk.interfaces.pl [-h|--help] [--no-if] list of file names\n";
    exit;
  }
  if($file =~ /^--no-preproc/) {
    $use_preproc=0;
    next;
  }
  if (! open(INPUT,"<$file") ) {
      print STDERR "Can't open input file $file.\n";
      next;}

  if($file =~ /^(.*)[.].*/) {
    print "#if ! defined(IN_$1)\n" if($use_preproc == 1);
    }

  while (<INPUT>) {
    if($_ =~ /^(.*)(!InTf!)(.*)/) {
      $line = $1 ;
      if($line  =~ /!!(.*)/) { 
        $line = $1;
      }   # lines starting with !!
      if($line  =~ /^!VOID[\$].*[ ]([a-zA-Z][a-zA-Z_0-9]*)[ ]*$/) {   # special IGNORE TKR
        if($use_preproc == 1) {
          print "# if __GNUC__  > 3 && __GNUC_MINOR__ > 7\n";           # gfortran 4.8 and up
          print" type(*), dimension(..) :: $1\n";                       # deferred type and rank
          print "#else\n";
          print" integer :: $1\n";
          print "# if __GNUC__  > 3 && __GNUC_MINOR__ > 8\n";           # gfortran 4.9 and up (never reached now)
          print"!GCC\$ ATTRIBUTES NO_ARG_CHECK :: $1\n";
          print "# elif defined(__INTEL_COMPILER)\n";                   # Intel ifort
          print"!DEC\$ ATTRIBUTES NO_ARG_CHECK :: $1\n";
          print "# elif defined (__PPC__)\n";                           # Power PC (IBM xlf)
          print"!IBM* ignore_tkr $1\n";
          print "# elif defined(__PGI)\n";                              # Portland Group Fortran
          print"!DIR\$ ignore_tkr $1\n";
          print "# else\n";
          print"!\$PRAGMA ignore_tkr $1\n";                             # Others (Sun Studio)
          print "# endif\n";
          print"#endif\n";
        }else{
          print" integer :: $1\n";
          print"!GCC\$ ATTRIBUTES NO_ARG_CHECK :: $1\n";
          print"!DEC\$ ATTRIBUTES NO_ARG_CHECK :: $1\n";
          print"!IBM* ignore_tkr $1\n";
          print"!DIR\$ ignore_tkr $1\n";
          print"!\$PRAGMA ignore_tkr $1\n";                             # Others (Sun Studio)
        }
        next;
      }   # !VOID$
      $line =~ s/^\s+// ; $line =~ s/\s+$// ;   # trim space at both ends
      $tail = '&' ; $head = "      " ;

      $l = length "$line";
      while ($l > 66) {                # fixed/free format split with continuation lines
        $body = substr $line, 0, 66;
        print "$head$body$tail\n" ;    # print first 66 characters of line with continuation
        $head = "     &";
        $line  = substr $line, 66 ;    # get rid of first 66 characters of line
        $l = length "$line";
      }
      print "$head$line\n";
    }   # !InTf!
  }   # <INPUT>

  print "#endif\n" if($use_preproc == 1);
}
