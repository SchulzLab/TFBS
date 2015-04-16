#! /bin/bash --

cut -f -4 $1 |
sed -e '/^\s*$/d' -e '/\s*#.*$/d'
