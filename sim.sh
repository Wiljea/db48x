#!/usr/bin/env bash
# DB48x simulator launcher — recomputes X DISPLAY, sets memory, 2x window
cd ~/db48x || exit 1
export DISPLAY=$(grep nameserver /etc/resolv.conf | awk "{print \$2}"):0
exec ./db48x -m 100000 -s2 "$@"
