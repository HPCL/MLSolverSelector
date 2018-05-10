#!/bin/bash --
../dependencies/sqlite/bin/sqlite3 -batch $1 <<"EOF"
.mode column
.headers on
SELECT * FROM DATA;
EOF
