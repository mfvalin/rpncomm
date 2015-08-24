#!/usr/bin/env perl

# Si pas d'arguments, on prend STDIN
if ( $#ARGV > -1 ){
 @listfile = (@ARGV) ;
 $listfile[$#listfile] =~ s/\n/ /;
}

foreach $file (@listfile) {

  if (! open(INPUT,"<$file") ) {
      print STDERR "Can't open input file $file.\n";
      next;}

  if($file =~ /^(.*)[.].*/) {
    print "#if defined(IN_$1)\n" ;
    }

  while (<INPUT>) {
    if($_ =~ /^(.*)(!InTf!)(.*)/) {
      $line = $1 ;
      if($line  =~ /!!(.*)/) { 
        $line = $1;
      }
      if($line  =~ /^!VOID[\$].*[ ]([a-zA-Z][a-zA-Z_0-9]*)[ ]*$/) {   # special IGNORE TKR
        print "# if __GNUC__  > 3 && __GNUC_MINOR__ > 7\n type(*), dimension(..) :: $1\n";
        print "#else\n integer :: $1\n";
        print "# if __GNUC__  > 3 && __GNUC_MINOR__ > 8\n!GCC\$ ATTRIBUTES NO_ARG_CHECK :: $1\n";
        print "# elif defined(__INTEL_COMPILER)\n!DEC\$ ATTRIBUTES NO_ARG_CHECK :: $1\n";
        print "# elif defined (__PPC__)\n!IBM* ignore_tkr $1\n";
        print "# elif defined(__PGI)\n!DIR\$ ignore_tkr $1\n";
        print "# else\n!\$PRAGMA ignore_tkr $1\n";
        print "# endif\n#endif\n";
        next;
      }
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
    }
  }

  print "#endif\n"
}
