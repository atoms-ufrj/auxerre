#!/bin/bash
echo Processing 1.pos
sed -E 's/(^[0-9]+)/Rb/g' 1.pos | sed -e '/^#/d' -e 's/Rb 10$/10\n/g' | gzip | pv > 1.xyz.gz
echo Processing 2.pos
sed -E 's/(^[0-9]+)/In/g' 2.pos | sed -e '/^#/d' -e 's/In 990$/990\n/g' | gzip | pv > 2.xyz.gz
