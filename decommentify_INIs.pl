#!/usr/bin/perl -w
use strict;

# decommentify_INI.pl
#
# A simple script to remove all comments (and blank lines) from an INI file.  This action can be undone with the script commentify_INI.pl.
#
# Josh Burton
# August 2013


# Check for input arguments.
if ( @ARGV != 1 ) {
    print STDERR "\nSyntax: $0 <INI-file>   (writes to stdout)\n\n";
    exit;
}


# Open the filename provided.
die "Can't find input file `$ARGV[0]`: $!\n" unless -e $ARGV[0];
open IN, '<', $ARGV[0] or die;

# Filter out empty lines and lines starting with a comment ('#').  Print all other lines.
map { print unless /^[\n\#]/ } <IN>;

close IN;
