#!/usr/bin/env bash
set -e

# This script starts a Docker container with Jupyter Lab running on port 8888 and prints the url to access it. It mounts the current working directory to the 'mnt' directory.
# It will delete a Docker container of the same name

port=8888
# container_tag="ont_simulator:v1"
container_tag="ghcr.io/ratschlab/sim_read_until:main_dev"
container_name="ont_simulator"

function usage_and_exit() {
    echo "Usage: start_docker_container.sh [-p port]"
    echo "Starts a Docker container with Jupyter Lab and prints the url to access it. It mounts the current working directory."
    echo "Options:"
    echo "-p port: The port to run Jupyter Lab on. Default: $port"
    echo "-t container_tag: The tag of the Docker container to run. Default: $container_tag"
    echo "-n container_name: The name of the Docker container to run. Default: $container_name"
    echo "-h: Print this help message and exit"
    exit 1
}

# make it possible to pass a port as an argument using getopt
while getopts "p:t:n:h" opt; do
  case $opt in
    p) port="$OPTARG"
    ;;
    t) container_tag="$OPTARG"
    ;;
    n) container_name="$OPTARG"
    ;;
    h) usage_and_exit
    ;;
    \?) echo "Invalid option -$opt" >&2
    usage_and_exit
    ;;
  esac
done

echo "Starting a Docker container with name '$container_name' (tag '$container_tag') on port $port with Jupyter Lab running, mounting the current working directory '$(pwd)'."

# container_name=$(docker run ...) # to get an unused name
docker rm -f "$container_name" || true
docker run --name "$container_name" -v "$(pwd):/root/ont_project_all/mnt" -p $port:$port -d -t "$container_tag" /bin/bash -c "cd /root/ont_project_all && source ~/ont_project_all/ont_project_venv/bin/activate && jupyter lab --port $port --allow-root --ip 0.0.0.0 --no-browser"
echo "Container name: $container_name"
sleep 3 # so that jupyter can startup and outputs the URL

# parse the url from the file
#     To access the server, open this file in a browser:
#         file:///root/.local/share/jupyter/runtime/jpserver-1-open.html
#     Or copy and paste one of these URLs:
#         http://ddb53c639c02:8888/lab?token=906627c0f1d9b6e1374c1a80095bf3a73e5190319232f48e
#         http://127.0.0.1:8888/lab?token=906627c0f1d9b6e1374c1a80095bf3a73e5190319232f48e
# [I 2023-07-18 14:19:36.095 ServerApp] Skipped non-installed server(s): bash-language-server, dockerfile-language-server-nodejs, javascript-typescript-langserver, jedi-language-server, julia-language-server, pyright, python-language-server, python-lsp-server, r-languageserver, sql-language-server, texlab, typescript-language-server, unified-language-server, vscode-css-languageserver-bin, vscode-html-languageserver-bin, vscode-json-languageserver-bin, yaml-language-server

# Given this output, parse the url "http://127.0.0.1:8888/lab?token=906627c0f1d9b6e1374c1a80095bf3a73e5190319232f48e" and open it in the browser.
output=$(docker logs "$container_name")

echo "Open the following URL to access Jupyter Lab:"
# get the line in the output matching "http://127.0.0.1"
echo "$output" | grep -F 'http://127.0.0.1'
# doing the following results in the URL being printed with the end concatenated, why? line_with_url=$(echo "$output" | grep -F 'http://127.0.0.1')

echo "If it does not contain the URL, inspect the output of 'docker logs $container_name' and find the URL."

echo >&2 "Once you are done, remove the container with 'docker rm -f $container_name'"