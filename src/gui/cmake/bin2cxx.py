# -*- coding: utf-8 -*-
import os
import re
import sys
from optparse import OptionParser

ROW_VALUES = 16

parser = OptionParser(usage="%prog [options] file")
parser.add_option("-i", "--identifier",
                  help="identifiers for the generated variables", metavar="IDENTIFIER")
parser.add_option("-o", "--output",
                  help="output filename", metavar="FILE")

(options, args) = parser.parse_args()

if len(args) != 1:
  sys.stderr.write ("Exactly one input file needs to be specified.\n")
  parser.print_help()
  sys.exit(1)
  
input_filename = args[0]
identifier = options.identifier
if not identifier:
  # No identifier arg, generate from input filename
  (input_dir, input_name) = os.path.split (input_filename)
  identifier = re.sub ("[^A-Za-z0-9_]", "_", input_name)
  if re.match("[0-9]", identifier): identifier = "_" + identifier
  
if options.output:
  output_file = open (options.output, "w")
else:
  output_file = sys.stdout

input_file = open (input_filename, "rb")
in_data = input_file.read()
in_len = len (in_data)
input_file.close()
  
output_file.write ("static const size_t %s_size = %d;\n" % (identifier, in_len))
output_file.write ("static const unsigned char %s_data[%s_size] = {" % (identifier, identifier))

n = 0
for byte in in_data:
  if (n % ROW_VALUES) == 0:
    output_file.write ("\n\t");
  else:
    output_file.write (" ")
  output_file.write ("0x%02x" % ord (byte))
  n = n+1
  if n < in_len: output_file.write (",")

output_file.write ("\n};\n")

