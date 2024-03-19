import os
import pandas as pd
import plotly.graph_objs as go

# Step 2: Load the Data
folder_data = '/Users/wfon4473/Documents/Bioinformatics/all_testdirs/testdir_heatmaps'


def clean_up(rpm_data):
    rpm_df = pd.read_csv(rpm_data, header=0, names=['taxon', "taxid", "rpm"])
    rpm_df = rpm_df[~rpm_df['taxon'].str.contains('No Hit')]
    return rpm_df

def merge_rpm_data(folder_data):
    all_sample_info = []
    for item in os.listdir(folder_data):
        if not item.startswith('.') and item.endswith('_rpm_data.csv'):
            sample_name = item.split("_rpm_data.csv")[0]
            sample_info = clean_up(os.path.join(folder_data, item))
            sample_info['sample_name'] = sample_name
            all_sample_info.append(sample_info)
    all_sample_info = pd.concat(all_sample_info).reset_index(drop=True)
    return all_sample_info

data = merge_rpm_data(folder_data)

# Step 3: Filter for Known Pathogens
# Load CSV containing known pathogens
known_pathogens = pd.read_csv('/Users/wfon4473/Documents/Bioinformatics/metagp_report/database/pathogen_list.csv')

# Filter data to include only known pathogens
pathogen_data = data[data['taxon'].isin(known_pathogens['Species']) | data['taxon'].isin(known_pathogens['AltNames'])]

def create_heatmap(df, title):
    heatmap_data = df.pivot(index='taxon', columns='sample_name', values='rpm')
    heatmap_data = heatmap_data.reindex(index=heatmap_data.sum(axis=1).sort_values(ascending=True).index)
    heatmap_height = len(heatmap_data) * 20
    heatmap_width = len(heatmap_data.columns) * 400
    heatmap = go.Figure(go.Heatmap(z=heatmap_data.values,
                                    x=heatmap_data.columns,
                                    y=heatmap_data.index,
                                    colorscale='OrRd'))
    heatmap.update_layout(title=title,
                          xaxis=dict(title='Sample', side='top', automargin=False, tickangle=-45),
                          yaxis=dict(title='Taxon', automargin=False),
                          height=heatmap_height,
                          width=1200,
                          plot_bgcolor='rgba(0,0,0,0)',
                          margin=dict(t=100, b=0, l=500, r=0),
                          xaxis_showgrid=False, yaxis_showgrid=False)
    return heatmap

heatmap_all = create_heatmap(data, 'Taxonomic Heatmap - All Taxa')
heatmap_pathogens = create_heatmap(pathogen_data, 'Taxonomic Heatmap - Known Pathogens')

# Export the HTML files
heatmap_all_html = heatmap_all.to_html(full_html=False)
heatmap_pathogens_html = heatmap_pathogens.to_html(full_html=False)

# Write the combined HTML content to a file
combined_html = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Combined Heatmaps</title>
    <script>
        function toggleHeatmap() {{
            var allHeatmap = document.getElementById('all-heatmap');
            var pathogensHeatmap = document.getElementById('pathogens-heatmap');
            var toggleButton = document.getElementById('toggle-button');
            if (allHeatmap.style.display === 'block') {{
                allHeatmap.style.display = 'none';
                pathogensHeatmap.style.display = 'block';
                toggleButton.textContent = 'Show All';
            }} else {{
                allHeatmap.style.display = 'block';
                pathogensHeatmap.style.display = 'none';
                toggleButton.textContent = 'Show Known Pathogens';
            }}
        }}
    </script>
</head>
<body>

<button id="toggle-button" onclick="toggleHeatmap()">Show Known Pathogens</button>

<div id="all-heatmap">{heatmap_all_html}</div>
<div id="pathogens-heatmap" style="display: none;">{heatmap_pathogens_html}</div>

</body>
</html>
"""

# Write the combined HTML content to a file
export_path = os.path.join(folder_data, "heatmap.html")
with open(export_path, "w") as file:
    file.write(combined_html)
