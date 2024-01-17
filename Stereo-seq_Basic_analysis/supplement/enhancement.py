import cv2
import numpy as np

##########################
# Histogram equalization #
##########################


def equalhist_transfer(
    img, method="global", cliplimit=20, tileGridSize=(7, 7)
) -> np.ndarray:
    """
    Histogram equalization for image enhancement, including:
        global histogram equalization,
        adaptive local histogram equalization.
    Args:
        img: Single-channel uint8 type grayscale image[0, 255].
        method: Histogram equalization methods. Available `method` are:
            * `'global'`: global histogram equalization;
            * `'local'`: adaptive local histogram equalization.
        cliplimit: Threshold to limit contrast when method = 'local'.
        tileGridSize: Size of grid for histogram equalization when method = 'local'.
    Returns:
        Image matrix after image enhancement.
    """
    if method == "global":
        # Global histogram equalization.
        return cv2.equalizeHist(img)
    elif method == "local":
        # Adaptive local histogram equalization.
        clahe = cv2.createCLAHE(clipLimit=cliplimit, tileGridSize=tileGridSize)
        return clahe.apply(img)
    else:
        raise ValueError(
            f"Available histogram equalization methods are `global` and `local`."
        )


####################
# Gamma correction #
####################


def gamma_correction(img, gamma=0.5):
    """
    Gamma correction is also known as the Power Law Transform.
    First, our image pixel intensities must be scaled from the range [0, 255] to [0, 1.0].
    From there, we obtain our output gamma corrected image by applying the following equation:
        Gamma values < 1 will shift the image towards the darker end of the spectrum;
        Gamma values > 1 will make the image appear lighter;
        Gamma value = 1 will have no affect on the input image.
    Args:
        img: Single-channel uint8 type grayscale image[0, 255].
        gamma: Gamma values.
               Gamma values < 1 will shift the image towards the darker end of the spectrum;
               Gamma values > 1 will make the image appear lighter.
               Gamma value = 1 will have no affect on the input image.
    Returns:
        Image matrix after image enhancement.
    Example:
        >>> gamma_correction(img, gamma=0.5)
    """
    img = 255 * np.power(img / 255, gamma)
    img = np.around(img)
    img[img > 255] = 255
    return img.astype(np.uint8)


#########################################
# Haze Removal Using Dark Channel Prior #
#########################################


def zmMinFilterGray(src, r=7):
    """Minimum Filtering."""
    return cv2.erode(src, np.ones((2 * r + 1, 2 * r + 1)))


def guidedfilter(I, p, r, eps):
    """Guided Filtering."""
    height, width = I.shape
    m_I = cv2.boxFilter(I, -1, (r, r))
    m_p = cv2.boxFilter(p, -1, (r, r))
    m_Ip = cv2.boxFilter(I * p, -1, (r, r))
    cov_Ip = m_Ip - m_I * m_p

    m_II = cv2.boxFilter(I * I, -1, (r, r))
    var_I = m_II - m_I * m_I

    a = cov_Ip / (var_I + eps)
    b = m_p - a * m_I

    m_a = cv2.boxFilter(a, -1, (r, r))
    m_b = cv2.boxFilter(b, -1, (r, r))
    return m_a * I + m_b


def Defog(src, min_r, guide_r, eps, w, maxV1):
    """Calculate the atmospheric mask image V1 and the light value A, V1 = 1-t/A"""
    V1 = np.min(src, 2)
    Dark_Channel = zmMinFilterGray(V1, min_r)

    V1 = guidedfilter(V1, Dark_Channel, guide_r, eps)
    bins = 2000
    ht = np.histogram(V1, bins)
    d = np.cumsum(ht[0]) / float(V1.size)
    for lmax in range(bins - 1, 0, -1):
        if d[lmax] <= 0.999:
            break
    A = np.mean(src, 2)[V1 >= ht[1][lmax]].max()
    V1 = np.minimum(V1 * w, maxV1)
    return V1, A


def deHaze(src, min_r=7, guide_r=81, eps=0.001, w=0.95, maxV1=0.80, bGamma=False):
    """
    References: Single Image Haze Removal Using Dark Channel Prior.
    Args:
        src: RGB image.
        min_r: Dark channel minimum filter radius. The larger the radius, the less obvious the effect of fog removal.
           The recommended range is generally between 5 and 25. Generally, 5, 7, and 9 will achieve good results.
        guide_r: The mean filter radius in guided filtering. The recommended value of this radius is not less than 4
                 times of the minimum filter radius when finding the dark channel. Because the dark channel after the
                 previous minimum value is block by block, in order to make the transmittance map more refined,
                 this r cannot be too small.
        eps: This value is only to ensure that the division sign is not 0, which is generally smaller.
             Basically no modification is required.
        w: The level of fog is preserved. Basically no modification is required.
        maxV1: Limit the range of image values.
        bGamma: Whether to use gamma correction at the end of image processing.
    """

    src_type = 255 if np.max(src) > 1 else 1
    if src_type == 255:
        src = src / 255

    Y = np.zeros(src.shape)
    Mask_img, A = Defog(src, min_r, guide_r, eps, w, maxV1)

    for k in range(3):
        Y[:, :, k] = (src[:, :, k] - Mask_img) / (1 - Mask_img / A)
    Y = np.clip(Y, 0, 1)
    if bGamma:
        Y = Y ** (np.log(0.5) / np.log(Y.mean()))

    if src_type == 255:
        Y = Y * 255
        Y[Y > 255] = 255

    return Y