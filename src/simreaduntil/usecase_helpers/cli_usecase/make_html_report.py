"""
Make an overview page from all figures in the directory
"""
import argparse
import logging
from pathlib import Path
import sys
from simreaduntil.shared_utils.logging_utils import add_comprehensive_stream_handler_to_logger, setup_logger_simple
from simreaduntil.shared_utils.utils import print_args

logger = setup_logger_simple(__name__)
"""module logger"""

mainpage_template = r"""
<!DOCTYPE html>
<html>
  <head>
    <title>Generated Figures</title>
  </head>
  <style>
    img {
      display: inline-block;
      margin: 5px;
      max-width: 100%;
      max-height: 400px;
      text-align: center;
    }
    
  </style>
  <body>
    <h1>Generated Figures</h1>
    
    {FIGURES_HTML}

  </body>
</html>
"""
figure_template = r"""
<h2 class="image-name">{figure_name}</h2>
    <div class="photo">
      <img src="{figure_filename}">
    </div>
    <hr style="height:2px;border-width:0;color:gray;background-color:gray">
    
"""

def parse_args(args=None):
    parser = argparse.ArgumentParser(description="Combine images into html report")
    parser.add_argument("figure_dir", type=Path, help="Directory containing figures (.png files)")
    
    args = parser.parse_args(args)
    print_args(args, logger=logger)
    return args

def main():
    log_level = logging.DEBUG
    logging.getLogger(__name__).setLevel(log_level) # log level of this script (running as __main__)
    add_comprehensive_stream_handler_to_logger(None, level=log_level)
    
    args = parse_args()
    
    figure_dir = args.figure_dir
    assert figure_dir.exists(), f"figure_dir '{figure_dir}' does not exist"
    
    figure_htmls = ""
    for filename in figure_dir.glob("*.png"):
        logger.info(f"Adding '{filename}'")
        # using relative path wrt to final file!
        figure_htmls += figure_template.format(figure_name=filename.name, figure_filename=filename.name)
    
    html_filename = figure_dir / "figures.html"
    with open(html_filename, "w") as f:
        f.write(mainpage_template.replace("{FIGURES_HTML}", figure_htmls))
    
  
# template
# <!DOCTYPE html>
# <html>
#   <head>
#     <title>Generated Figures</title>
#   </head>
#   <style>
#     img {
#       display: inline-block;
#       margin: 5px;
#       width: 100%;
#       max-height: 400px;
#       text-align: center;
#     }
    
#   </style>
#   <body>
#     <h1>Generated Figures</h1>
    
#     <h2 class="image-name">nb_basepairs_recrej_per_channel.png</h2>
#     <div class="photo">
#       <img src="figures/nb_basepairs_recrej_per_channel.png">
#     </div>
#     <hr style="height:2px;border-width:0;color:gray;background-color:gray">

#   </body>
# </html>