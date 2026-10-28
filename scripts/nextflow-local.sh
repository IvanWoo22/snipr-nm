#!/usr/bin/env bash
# Offline launcher for this project's existing Nextflow 25.10.4 / Java 21 runtime.
# Module-open options follow the official v25.10.4 launcher (lines 423-439):
# https://github.com/nextflow-io/nextflow/blob/v25.10.4/nextflow
# Running java -jar without these options prevents cache serialization on Java 21.
set -euo pipefail
SNIPR_PROJECT=$(cd "$(dirname "$(realpath "$0")")/.." && pwd)
SNIPR_NXF_HOME=${NXF_HOME:-"$SNIPR_PROJECT/.runtime/nextflow"}
SNIPR_JAVA=${JAVA_HOME:+"$JAVA_HOME/bin/"}java
SNIPR_JAR="$SNIPR_NXF_HOME/framework/25.10.4/nextflow-25.10.4-one.jar"
if [[ ! -f "$SNIPR_JAR" ]]; then
    echo "Missing Nextflow runtime: $SNIPR_JAR. Install the official Nextflow distribution." >&2
    exit 1
fi
SNIPR_JAVA_FLAGS=()
for module in java.lang java.io java.nio java.net java.util java.util.concurrent.locks \
    java.util.concurrent.atomic java.nio.file.spi sun.nio.ch sun.nio.fs \
    sun.net.www.protocol.http sun.net.www.protocol.https sun.net.www.protocol.ftp \
    sun.net.www.protocol.file jdk.internal.misc jdk.internal.vm java.util.regex; do
    SNIPR_JAVA_FLAGS+=("--add-opens=java.base/$module=ALL-UNNAMED")
done
export NXF_HOME="$SNIPR_NXF_HOME"
exec "$SNIPR_JAVA" "${SNIPR_JAVA_FLAGS[@]}" -jar "$SNIPR_JAR" "$@"
