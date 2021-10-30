#!/usr/bin/awk -f

/href/ { split($0, s, "\""); print "http://www.netlib.org/port/"s[2] }
