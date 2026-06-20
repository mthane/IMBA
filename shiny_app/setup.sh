#!/bin/bash
# Run from repository root:  bash shiny_app/setup.sh
# Or:  cd shiny_app && bash setup.sh
set -euo pipefail
cd "$(dirname "$0")"

sudo apt install build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libcairo2-dev
sudo apt install cmake

sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/'
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo apt update && sudo apt upgrade
sudo apt install r-base

echo "Package root: $PWD"
echo "Start R and run:"
echo "  install.packages('devtools')"
echo "  devtools::install()"
echo "  library(imba); imba::run_app()"
sudo -i R
