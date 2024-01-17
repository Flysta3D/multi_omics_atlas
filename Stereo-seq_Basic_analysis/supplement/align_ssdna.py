from typing import Optional, Tuple

import cv2
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.path import Path
from matplotlib.widgets import LassoSelector

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal


def read_bgi_as_dataframe(path: str) -> pd.DataFrame:
    """Read a BGI read file as a pandas DataFrame.

    Args:
        path: Path to read file.

    Returns:
        Pandas Dataframe with column names `gene`, `x`, `y`, `total` and
        additionally `spliced` and `unspliced` if splicing counts are present.
    """
    return pd.read_csv(
        path,
        sep="\t",
        encoding = 'unicode_escape',
        dtype={
            "geneID": "category",  # geneID
            "x": np.uint32,  # x
            "y": np.uint32,  # y
            "MIDCounts": np.uint16,  # total
        },
        comment="#",
    )


def equalhist_transfer(img, method="global", cliplimit=20) -> np.ndarray:
    """
    Histogram equalization for image enhancement, including:
        global histogram equalization,
        adaptive local histogram equalization.

    Args:
        img: Single-channel uint8 type grayscale image[0, 255].
        method: Histogram equalization methods. Available `method` are:
            * `'global'`: global histogram equalization;
            * `'local'`: adaptive local histogram equalization.
        cliplimit: Threshold to limit contrast when method = 'local' .
    Returns:
        Image matrix after image enhancement.
    """
    if method == "global":
        # Global histogram equalization.
        return cv2.equalizeHist(img)
    elif method == "local":
        # Adaptive local histogram equalization.
        clahe = cv2.createCLAHE(clipLimit=cliplimit, tileGridSize=(7, 7))
        return clahe.apply(img)
    else:
        raise ValueError(
            f"Available histogram equalization methods are `global` and `local`."
        )


def output_img(
    img: np.ndarray,
    filename: str = None,
    window_size: tuple = (1024, 1024),
    show_img: bool = True,
):
    """
    Output the image matrix as a 2D image file.

    Args:
        img: Image matrix.
        filename: Output image filename, the end of which can be .bmp, .dib, .jpeg, .jpg, .jpe, .png, .webp, .pbm,
                  .pgm, .ppm, .pxm, .pnm, .sr, .ras, .tiff, .tif, .exr, .hdr, .pic, etc.
        window_size: The size of the image visualization window.
        show_img: Whether to create a window to display the image.
    """
    if show_img:
        cv2.namedWindow("Image", cv2.WINDOW_NORMAL | cv2.WINDOW_KEEPRATIO)
        cv2.resizeWindow("Image", window_size[0], window_size[1])
        cv2.imshow("Image", img)
    if filename is not None:
        cv2.imwrite(filename=filename, img=img)
'''    
    cv2.waitKey(0)
    cv2.destroyAllWindows()

    if show_img:
        cv2.namedWindow("Image", cv2.WINDOW_NORMAL | cv2.WINDOW_KEEPRATIO)
        cv2.resizeWindow("Image", window_size[0], window_size[1])
        cv2.imshow("Image", img)
    if filename is not None:
        cv2.imwrite(filename=filename, img=img)
    cv2.waitKey(0)
    cv2.destroyWindow("Image")
    cv2.waitKey(1)
'''

class SelectFromCollection:
    """
    Select indices from a matplotlib collection using `LassoSelector`.

    Selected indices are saved in the `ind` attribute. This tool fades out the
    points that are not part of the selection (i.e., reduces their alpha
    values). If your collection has alpha < 1, this tool will permanently
    alter the alpha values.

    Note that this tool selects collection objects based on their *origins*
    (i.e., `offsets`).

    Parameters
    ----------
    ax : `~matplotlib.axes.Axes`
        Axes to interact with.
    collection : `matplotlib.collections.Collection` subclass
        Collection you want to select from.
    alpha_other : 0 <= float <= 1
        To highlight a selection, this tool sets all selected points to an
        alpha value of 1 and non-selected points to *alpha_other*.
    """

    def __init__(self, ax, collection, alpha_other=0.3):
        self.canvas = ax.figure.canvas
        self.collection = collection
        self.alpha_other = alpha_other

        self.xys = collection.get_offsets()
        self.Npts = len(self.xys)

        # Ensure that we have separate colors for each object
        self.fc = collection.get_facecolors()
        if len(self.fc) == 0:
            raise ValueError("Collection must have a facecolor")
        elif len(self.fc) == 1:
            self.fc = np.tile(self.fc, (self.Npts, 1))

        self.lasso = LassoSelector(
            ax, onselect=self.onselect, props=dict(color="black")
        )
        self.ind = []

    def onselect(self, verts):
        path = Path(verts)
        self.ind = np.nonzero(path.contains_points(self.xys))[0]
        self.fc[:, -1] = self.alpha_other
        self.fc[self.ind, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()

    def disconnect(self):
        self.lasso.disconnect_events()
        self.fc[:, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()


def lasso_2d(
    coords: np.ndarray,
    values: Optional[np.ndarray] = None,
    cmap: str = "viridis",
    point_size: int = 10,
) -> np.ndarray:
    """
    Pick the interested part of coordinates by interactive approach.

    Args:
        coords: The coordinates of points data.
                [[x1, y1],
                 [x2, y2],
                 [x3, y3],
                 [x4, y4]
                 ...]
        values: The values of points data.
                [value1, value2, value3, value4,...]
        cmap: Visualization colors for point data.
        point_size: Points size.
    Returns:
        select_points: The coordinates of the last selected points.
    """

    def fig_set():
        mpl.use("TkAgg")
        fig1, ax1 = plt.subplots(
            subplot_kw=dict(
                xlim=(coords[:, 0].min(), coords[:, 0].max()),
                ylim=(coords[:, 1].min(), coords[:, 1].max()),
                autoscale_on=False,
            )
        )
        pts1 = ax1.scatter(
            coords[:, 0], coords[:, 1], s=point_size, c=values, cmap=cmap
        )
        ax1.set_title("Press enter to accept selected points.")
        return fig1, ax1, pts1

    def select(event):
        if event.key == "enter":
            selector.disconnect()
            ax.set_title("The points have been selected, please exit the window.")
            fig.canvas.draw()

    # Create lasso window.
    fig, ax, pts = fig_set()
    selector = SelectFromCollection(ax, pts)
    fig.canvas.mpl_connect("key_press_event", select)
    plt.show()

    # Output the coordinates of the selected points.
    select_points = selector.xys[selector.ind].astype(coords.dtype)
    return select_points


def lasso_pre(
    path: str,
    lasso: bool = True,
    EH_method: Literal["global", "local"] = "global",
    EH_cliplimit: int = 20,
    color_flip: bool = False,
    gray_factor: int = 1,
    show: bool = True,
    save_img: Optional[str] = None,
    save_lasso: Optional[str] = None,
) -> Tuple[pd.DataFrame, np.ndarray]:
    """
    Preprocessing of original lasso data.

    Args:
        path: Path to read file.
        lasso: Whether to lasso through the interactive window.
        EH_method: Histogram equalization methods. Available `method` are:
            * `'global'`: global histogram equalization;
            * `'local'`: adaptive local histogram equalization.
        EH_cliplimit: Threshold to limit contrast when method = 'local' .
        color_flip: If color_flip is True, flip the color of the image.
        gray_factor:  Increasing the value in the grayscale image.
        show: If show is True, generate a visual window to display the image.
        save_img: If save is not None, save is the path to save the image.
        save_lasso: If save is not None, save is the path to save the lasso data.
    """

    # Read bgi file.
    data = read_bgi_as_dataframe(path=path)

    # Lasso bgi data.
    coords_data = (
        data[["x", "y", "MIDCounts"]]
        .groupby(["x", "y"])["MIDCounts"]
        .sum()
        .to_frame("MIDCounts")
        .reset_index()
    )

    if lasso:
        select_coords = lasso_2d(
            coords=coords_data[["x", "y"]].values,
            values=coords_data["MIDCounts"].values,
        )
        select_coords = pd.DataFrame(select_coords, columns=["x", "y"])
        lasso_data = pd.merge(data, select_coords, on=["x", "y"], how="inner")
    else:
        select_coords = coords_data.copy()
        lasso_data = data.copy()

    if save_lasso is not None:
        lasso_data.to_csv(save_lasso, sep="\t", index=False)

    # Create lasso image matrix.
    img = pd.pivot_table(
        select_coords, index=["x"], columns=["y"], values="MIDCounts", fill_value=0
    ).values.astype(np.uint8)

    #  Increasing the value in the grayscale image.
    img = img * gray_factor
    img[img > 255] = 255

    # Image enhancement based on Histogram equalization.
    img = equalhist_transfer(img, method=EH_method, cliplimit=EH_cliplimit)

    # Flip the color of the image.
    img = cv2.bitwise_not(src=img) if color_flip else img

    output_img(img=img, filename=save_img, show_img=show)

    return lasso_data, img


def crop_by_ssdna(
    path, path_ssdna, save_lasso=None, save_img=None, show=True
) -> Tuple[pd.DataFrame, np.ndarray]:
    """
    Cropping based on ssDNA.
    Args:
        path: Path to read lasso file.
        path_ssdna: Path to read ssDNA image file (grayscale).
        save_lasso: If save is not None, save is the path to save the lasso data.
        save_img: If save is not None, save is the path to save the image.
        show: If show is True, generate a visual window to display the image.
    """

    #data = read_bgi_as_dataframe(path=path)
    data = pd.read_csv(path,comment="#",sep="\t",encoding = 'unicode_escape')
    
    img = cv2.imread(path_ssdna, 0)
    background = img[0, 0]

    x_list = np.arange(start=data["x"].min(), stop=data["x"].max() + 1, step=1)
    y_list = np.arange(start=data["y"].min(), stop=data["y"].max() + 1, step=1)

    img_table = pd.DataFrame(img, index=x_list, columns=y_list)
    img_table["x"] = img_table.index
    img_data = pd.melt(img_table, id_vars=["x"])
    img_data.columns = ["x", "y", "value"]
    img_data = img_data[img_data["value"] != background]

    cropped_data = pd.merge(data, img_data[["x", "y"]], on=["x", "y"], how="inner")
    if save_lasso is not None:
        cropped_data.to_csv(save_lasso, sep="\t", index=False)

    cropped_x_list = np.arange(
        start=cropped_data["x"].min(), stop=cropped_data["x"].max() + 1, step=1
    )
    cropped_y_list = np.arange(
        start=cropped_data["y"].min(), stop=cropped_data["y"].max() + 1, step=1
    )
    cropped_img = img_table.loc[cropped_x_list, cropped_y_list]
    if save_img is not None:
        output_img(img=np.asarray(cropped_img), filename=save_img, show_img=show)

    return cropped_data, cropped_img
