import os
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score


class BioRefineryModel:

    def __init__(self, folder_path, product_name, sizes = None,
                 sheet_name='MSP_state_all',
                 USA_RAINFED=None, GREAT_LAKES=None):

        self.folder_path = folder_path
        self.product_name = product_name
        self.sizes = sizes
        self.sheet_name = sheet_name
        self.output_folder = os.path.join(self.folder_path, "outputs")
        os.makedirs(self.output_folder, exist_ok=True)
        # If not provided, import from project.paths
        if USA_RAINFED is None or GREAT_LAKES is None:
            try:
                from project.paths import US_RAINFED as USA_RF, GREAT_LAKES as GL
                self.USA_RAINFED = USA_RF
                self.GREAT_LAKES = GL
            except ImportError:
                raise ImportError(
                    "USA_RAINFED and GREAT_LAKES not provided and couldn't import from project.paths"
                )
        else:
            self.USA_RAINFED = USA_RAINFED
            self.GREAT_LAKES = GREAT_LAKES

    # =========================
    # 1. LOAD DATA
    # =========================
    def _process_files(self, tag):

        sizes, mean_list, median_list = [], [], []

        for filename in os.listdir(self.folder_path):
            if filename.endswith('.xlsx') and (tag in filename):

                filepath = os.path.join(self.folder_path, filename)

                size = filename.split('_')[0]
                try:
                    size = float(size)
                except ValueError:
                    continue

                df = pd.read_excel(filepath, sheet_name=self.sheet_name, index_col=0)
                df = df.apply(pd.to_numeric, errors='coerce')

                sizes.append(size)
                mean_list.append(df.mean(axis=0).values)
                median_list.append(df.median(axis=0).values)

        mean_df = pd.DataFrame(mean_list, columns=df.columns)
        median_df = pd.DataFrame(median_list, columns=df.columns)

        mean_df.insert(0, 'size', sizes)
        median_df.insert(0, 'size', sizes)

        mean_df = mean_df.sort_values(by='size').reset_index(drop=True)
        median_df = median_df.sort_values(by='size').reset_index(drop=True)

        # Replace with real sizes
        # Only override sizes if provided
        if self.sizes is not None:
            if len(self.sizes) != len(mean_df):
                raise ValueError("Provided sizes do not match number of files")

            mean_df['size'] = self.sizes
            median_df['size'] = self.sizes

        return mean_df, median_df

    # =========================
    # 2. A COEFFICIENT
    # =========================
    def _compute_a_coeff(self, df_p0, df_p1, label):

        assert all(df_p0['size'] == df_p1['size'])
        assert all(df_p0.columns == df_p1.columns)

        values = df_p1.drop(columns='size').values - df_p0.drop(columns='size').values

        a_df = pd.concat(
            [df_p0[['size']],
             pd.DataFrame(values, columns=df_p0.columns[1:])],
            axis=1
        )

        stats = a_df.drop(columns='size')

        print(f"\n--- {label} a-coefficient stats ({self.product_name}) ---")
        print("Min:", stats.min().min())
        print("Max:", stats.max().max())
        print("Range:", stats.max().max() - stats.min().min())
        print("Std:", stats.stack().std())

        return a_df

    # =========================
    # 3. FIT MODELS
    # =========================
    def _fit_models(self, df):

        def sat_2(x, a, b):
            return a/x + b

        models = {"saturation-2": (sat_2, ["a", "b"])}

        x = df['size'].values
        results = []

        for state in df.columns[1:]:
            y = df[state].values

            for name, (func, params) in models.items():
                try:
                    popt, _ = curve_fit(func, x, y, maxfev=10000)
                    y_pred = func(x, *popt)
                    r2 = r2_score(y, y_pred)

                    result = {"NAME": state, "Model": name, "R2": r2}
                    for i, p in enumerate(params):
                        result[p] = popt[i]

                    results.append(result)

                except Exception:
                    continue

        results_df = pd.DataFrame(results).round(3)
        results_df = results_df.sort_values(by=["NAME", "R2"], ascending=[True, False])

        # output_name = f"results_{self.product_name}_df_mean.csv"
        # results_df.to_csv(output_name, index=False)
        output_path = os.path.join(self.output_folder, f"results_{self.product_name}_df_mean.csv")
        results_df.to_csv(output_path, index=False)
        print(f"\nSaved fit results file to: {output_path}")

        return results_df

    # =========================
    # 4. SPATIAL SAMPLING
    # =========================
    def _generate_spatial_data(self, results_df):

        USA = gpd.read_file(self.USA_RAINFED).to_crs('EPSG:5070')
        GL = gpd.read_file(self.GREAT_LAKES).to_crs('EPSG:5070')

        USA = gpd.overlay(USA, GL, how="difference")

        USA['NAME'] = USA['NAME'].str.upper()
        results_df['NAME'] = results_df['NAME'].str.upper()

        merged = USA.merge(results_df, on="NAME")

        xmin, ymin, xmax, ymax = USA.total_bounds

        np.random.seed(123)
        N = 1_000_000

        lats = np.random.uniform(ymin, ymax, N)
        lons = np.random.uniform(xmin, xmax, N)

        points = gpd.GeoDataFrame(
            geometry=[Point(lon, lat) for lon, lat in zip(lons, lats)],
            crs=USA.crs
        )

        points = gpd.sjoin(points, merged[["NAME", "geometry", "a", "b"]],
                           how="left", predicate="within")

        points["a"] = points["a"].fillna(results_df.a.max())
        points["b"] = points["b"].fillna(results_df.b.max())

        # Save outputs
        # points.to_csv(f"state_samples_{self.product_name}.csv", index=False)
        csv_path = os.path.join(self.output_folder, f"state_samples_{self.product_name}.csv")
        points.to_csv(csv_path, index=False)
        print(f"Saved state samples to: {csv_path}")

        X_data = np.column_stack((lons, lats))
        # np.save(f"X_data_{self.product_name}.npy", X_data)
        x_path = os.path.join(self.output_folder, f"X_data_{self.product_name}.npy")
        np.save(x_path, X_data)
        print(f"Saved x data train to: {x_path}")

        outputs = np.column_stack((points["a"], points["b"]))
        # np.save(f"outputs_{self.product_name}.npy", outputs)
        outputs_path = os.path.join(self.output_folder, f"outputs_{self.product_name}.npy")
        np.save(outputs_path, outputs)
        print(f"Saved b and c coefficients to: {outputs_path}")

    # =========================
    # MAIN PIPELINE
    # =========================
    def run_all(self):
        print(f"\n=== Running analysis for {self.product_name} ===")

        mean_p0, median_p0 = self._process_files('_p0')
        mean_p1, median_p1 = self._process_files('_p1')

        self._compute_a_coeff(mean_p0, mean_p1, "Mean")
        self._compute_a_coeff(median_p0, median_p1, "Median")

        results_df = self._fit_models(mean_p0)

        if self.USA_RAINFED and self.GREAT_LAKES:
            self._generate_spatial_data(results_df)
            
        print(f"\n=== Finished analysis for {self.product_name} ===")
            
## Example usage:
# analysis = BioRefineryModel(
#     folder_path=folder_path,
#     product_name="ethanol",
#     sizes=[56.5, 80.7, 151.3, 252.05, 504.05]
# )

# analysis.run_all()