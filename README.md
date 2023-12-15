# SpectrumViewer
Code for viewing data cubes from the ALFALFA survey. The datacubes were going up on the NRAO's website so it was suggested that I make an interactive viewer for the datacubes. The user would input a data cube, and will presented with a window with two panels, where on the left panel users have an image of the sky(RA on x axis and DEC on y) in radio (extracted from the data cube) and a slider that lets them scroll through different frequencies and updates the image in real time. And when users spot a galaxy or anything interesting, they can click on it and it generates a spectrum with flux on y axis and frequency on x axis so the spectrum could be analyzed to see if there was an HI profile or other interesting features and then if the user wants they can save the spectrum as a FITS file.

# Dependencies
- numpy
- matplotlib
- astropy
- matplotlib.widgets
- matplotlib.gridspec

To install the required dependencies, use pip:

''' bash
pip install numpy matplotlib astropy'''

# Usage
from spectrum_viewer import SpectrumViewer

*Initialize the viewer with a FITS file*

viewer = SpectrumViewer(fits_file='your_file.fits')

*Use the main method to display an inital slice*
viewer.main(channel)

For more specifics, see SpecViewer.py. I have commented it extensively so it shouldn't be too hard to get a grasp on various functionalities.

