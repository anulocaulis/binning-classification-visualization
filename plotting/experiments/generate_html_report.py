#!/usr/bin/env python3
"""Generate HTML report of all experimental sandbox figures."""

import os
import sys
import base64
from pathlib import Path
from datetime import datetime

def generate_html_report(sandbox_dir: str = "plots/sandbox", outpath: str = "experimental_figures.html"):
    """Generate HTML report from PNG files in sandbox directory.
    
    Args:
        sandbox_dir: Path to sandbox directory containing PNG files
        outpath: Output HTML file path
    """
    sandbox_path = Path(sandbox_dir)
    
    if not sandbox_path.exists():
        print(f"ERROR: Sandbox directory not found: {sandbox_dir}")
        return
    
    # Collect all PNG files
    png_files = sorted(sandbox_path.glob("*.png"))
    
    if not png_files:
        print(f"WARNING: No PNG files found in {sandbox_dir}")
        return
    
    # Categorize files
    categories = {
        "Contig Counts": [],
        "Rarefaction Curves": [],
        "Ng50 / L50 Metrics": [],
        "Other": [],
    }
    
    for png_file in png_files:
        name = png_file.name
        if "contig_counts_triptych" in name:
            categories["Contig Counts"].append(png_file)
        elif "rarefaction" in name:
            categories["Rarefaction Curves"].append(png_file)
        elif "ng50" in name or "l50" in name:
            categories["Ng50 / L50 Metrics"].append(png_file)
        else:
            categories["Other"].append(png_file)
    
    # Generate HTML
    html_lines = [
        "<!DOCTYPE html>",
        "<html lang='en'>",
        "<head>",
        "  <meta charset='UTF-8'>",
        "  <meta name='viewport' content='width=device-width, initial-scale=1.0'>",
        "  <title>Experimental Figures - Binning & Classification Visualization</title>",
        "  <style>",
        "    body {",
        "      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;",
        "      line-height: 1.6;",
        "      color: #333;",
        "      background: #f5f5f5;",
        "      margin: 0;",
        "      padding: 20px;",
        "    }",
        "    .container {",
        "      max-width: 1400px;",
        "      margin: 0 auto;",
        "      background: white;",
        "      padding: 40px;",
        "      border-radius: 8px;",
        "      box-shadow: 0 2px 4px rgba(0,0,0,0.1);",
        "    }",
        "    h1 {",
        "      color: #2c3e50;",
        "      border-bottom: 3px solid #3498db;",
        "      padding-bottom: 15px;",
        "      margin-bottom: 10px;",
        "    }",
        "    .metadata {",
        "      font-size: 0.9em;",
        "      color: #7f8c8d;",
        "      margin-bottom: 30px;",
        "    }",
        "    h2 {",
        "      color: #34495e;",
        "      margin-top: 40px;",
        "      margin-bottom: 20px;",
        "      border-left: 4px solid #3498db;",
        "      padding-left: 15px;",
        "    }",
        "    .figure-grid {",
        "      display: grid;",
        "      grid-template-columns: repeat(auto-fit, minmax(500px, 1fr));",
        "      gap: 30px;",
        "      margin-bottom: 40px;",
        "    }",
        "    .figure {",
        "      border: 1px solid #ecf0f1;",
        "      border-radius: 6px;",
        "      overflow: hidden;",
        "      background: #fafafa;",
        "      box-shadow: 0 1px 3px rgba(0,0,0,0.08);",
        "      transition: box-shadow 0.3s ease;",
        "    }",
        "    .figure:hover {",
        "      box-shadow: 0 4px 8px rgba(0,0,0,0.12);",
        "    }",
        "    .figure img {",
        "      width: 100%;",
        "      height: auto;",
        "      display: block;",
        "    }",
        "    .figure-caption {",
        "      padding: 12px 15px;",
        "      font-size: 0.85em;",
        "      color: #555;",
        "      background: #f9f9f9;",
        "      border-top: 1px solid #ecf0f1;",
        "    }",
        "    .figure-name {",
        "      font-weight: 600;",
        "      color: #2c3e50;",
        "    }",
        "    .empty-section {",
        "      color: #95a5a6;",
        "      font-style: italic;",
        "      padding: 20px;",
        "      background: #ecf0f1;",
        "      border-radius: 4px;",
        "    }",
        "    .stats {",
        "      display: flex;",
        "      gap: 30px;",
        "      margin-top: 30px;",
        "      padding-top: 30px;",
        "      border-top: 2px solid #ecf0f1;",
        "      font-size: 0.9em;",
        "    }",
        "    .stat-item {",
        "      flex: 1;",
        "    }",
        "    .stat-label {",
        "      color: #7f8c8d;",
        "      font-size: 0.85em;",
        "      text-transform: uppercase;",
        "      letter-spacing: 0.5px;",
        "    }",
        "    .stat-value {",
        "      font-size: 1.3em;",
        "      font-weight: 600;",
        "      color: #3498db;",
        "      margin-top: 5px;",
        "    }",
        "  </style>",
        "</head>",
        "<body>",
        "  <div class='container'>",
        f"    <h1>Experimental Figures</h1>",
        f"    <p class='metadata'>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>",
        f"    <p class='metadata'>Source: {sandbox_path.resolve()}</p>",
    ]
    
    # Add each category
    total_figures = 0
    for category, files in categories.items():
        if not files:
            continue
        
        total_figures += len(files)
        html_lines.append(f"    <h2>{category}</h2>")
        html_lines.append("    <div class='figure-grid'>")
        
        for png_file in files:
            # Encode image as base64
            with open(png_file, 'rb') as img_f:
                img_data = base64.b64encode(img_f.read()).decode('utf-8')
            
            html_lines.append("      <div class='figure'>")
            html_lines.append(f"        <img src='data:image/png;base64,{img_data}' alt='{png_file.name}'>")
            html_lines.append("        <div class='figure-caption'>")
            html_lines.append(f"          <div class='figure-name'>{png_file.stem.replace('_', ' ').title()}</div>")
            html_lines.append(f"          <div>{png_file.name}</div>")
            html_lines.append("        </div>")
            html_lines.append("      </div>")
        
        html_lines.append("    </div>")
    
    # Add statistics
    html_lines.append("    <div class='stats'>")
    html_lines.append("      <div class='stat-item'>")
    html_lines.append("        <div class='stat-label'>Total Figures</div>")
    html_lines.append(f"        <div class='stat-value'>{total_figures}</div>")
    html_lines.append("      </div>")
    html_lines.append("      <div class='stat-item'>")
    html_lines.append("        <div class='stat-label'>Categories</div>")
    html_lines.append(f"        <div class='stat-value'>{sum(1 for files in categories.values() if files)}</div>")
    html_lines.append("      </div>")
    html_lines.append("    </div>")
    
    html_lines.extend([
        "  </div>",
        "</body>",
        "</html>",
    ])
    
    # Write HTML file
    html_content = "\n".join(html_lines)
    with open(outpath, "w") as f:
        f.write(html_content)
    
    print(f"✓ Generated: {outpath}")
    print(f"  - {total_figures} figures across {sum(1 for files in categories.values() if files)} categories")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sandbox-dir", default="plots/sandbox", help="Sandbox directory with PNG files")
    parser.add_argument("--outpath", default="experimental_figures.html", help="Output HTML file")
    args = parser.parse_args()
    
    generate_html_report(sandbox_dir=args.sandbox_dir, outpath=args.outpath)
