#! /bin/bash --

cut -f -3 $1 |
sed -e '/^\s*$/d' -e '/\s*#.*$/d'
