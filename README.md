# **deepNGS Navigator Panel Application**

## **Getting Started**

### 1. **Environment Setup**
Set up the required environment by creating a new Conda environment using the provided `environment.yml` file:

```bash
conda env create -f environment.yml
```

### 2. **Input Dataset**
The application requires the output of deepNGS mappings with the following minimum required columns:

- **`e1` and `e2`**: Two-dimensional embeddings used for visualization.
- **`AA`**: Sequence information that will be used by the panel application to draw MSAs.
- **`picked_clones`**: A string identifying clones of interest (e.g., binders). If no such clones are identified, this column can be left null.

#### **Additional Columns**:
Other columns in the dataset can be used for point size and color annotations.

#### **Setting Up Input Files**:
To integrate input files into the panel application, include their paths and metadata in the `processed_files.csv` table. This file should contain the following columns:
- **`name`**: A string identifier for the project. Use the format `XXX:YYY` if you want several subprojects (e.g., `YYY`) grouped under a main project (e.g., `XXX`) in the panel menu.
- **`path`**: The file path to the dataset corresponding to the project or subproject.

### 3. **Run the Panel Application**
Start the panel application using the following command:

```bash
panel serve main.py --port 5018
```

### 4. **Using the Application**
Once the panel application is running:
1. Copy the HTTP address provided in the terminal into your web browser.
2. In the browser interface:
   - Select the desired project from the menu.
   - Choose options for point size and color.
   - Press the **Display Map** button.

3. After the map is displayed:
   - Use the zoom controls to explore different parts of the map.
   - Use lasso or box select tools to select sequences of interest.
   - View their MSA patterns.
   - If needed, download the selected sequences directly from the interface.

Enjoy exploring your deepNGS data with the deepNGS Navigator Panel application!