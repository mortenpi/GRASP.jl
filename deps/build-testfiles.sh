#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd ${DIR}/../test/grasp/mixing || exit 1

$GRASP/bin/rangular <<-EOF
	y
EOF

$GRASP/bin/rwfnestimate <<-EOF
	y
	2
	*
EOF

$GRASP/bin/rmcdhf <<-EOF
	y
	1-3
	1
	1-2
	5
	*
	*
	100
EOF
