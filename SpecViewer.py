import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import matplotlib.widgets as widgets
import matplotlib.gridspec as gridspec
from astropy.constants import c
from astropy.table import Table
class SpectrumViewer:
    def __init__(self, filename= 'f',fits_file = '1220+09a_spectral.fits'):
        """
        Initialize the SpectrumViewer with a FITS file. 
        
        Parameters:
            fits_file (str): Path to the FITS file containing the data cube.

        This just means whenever you call this 'Object' or 'Class', which is just the big container for all the data processing and plotting, 
        it creates some variables and other stuff. It sets that to some arbitrary initial values. We do this so we can use the same variable
        across different methods or functions.
       
         As an example, say I have a Class called X, and in there I initialize two variables a = 0 and b = 0. Then I have two methods(functions),
        Method 1 and Method 2, I can have Method 1 say a = 1 and then in Method 2 I can say return a + b and it will return 1 + 0 = 1.
        Saving me the trouble of defining a = 1 again in Method 2. This is a very simple example but it's the same idea.
        """
        # Load the FITS file and store the data cube and header
        # with fits.open(fits_file) as hdul: is a way to open the file and close it automatically when you're done with it. 
        # It is the same as doing hdul = fits.open(fits_file) and then hdul.close() at the end of the code. 
        with fits.open(fits_file) as hdul:
            self.data_cube = hdul[0].data
            self.header = hdul[0].header

            
        # Create a unified grid layout for the figure. GridSpec makes a grid of 1 row and 2 columns. This is a blank canvas we can add stuff to.
        # Allocate 3 parts to the slice and 3 part to the spectrum.
        self.fig = plt.figure(figsize=(12, 8))
        gs = gridspec.GridSpec(1, 2, width_ratios=[5, 5])
        # Create axes for the slice and the spectrum in this unified layout.
        self.slice_ax = plt.subplot(gs[0])
        self.spectrum_ax = plt.subplot(gs[1])
        self.slice_ax.get_xlim()
        # Initialize the spectrum plot with an empty line.
        self.spectrum_line, = self.spectrum_ax.plot([], [])
        
        # Average the data cube along the polarization axis to get a 3D data cube. cube = [channel, y, x] where each combination of the three has a flux value.
        self.avgcube = np.mean(self.data_cube, axis=0)
        
        #Doing this so when user clicks, we just update these initialized figures rather than creating new figures for every click.
        self.spectrum_line, = self.spectrum_ax.plot([], [])
        # Create a slider axis at the bottom of the slice figure
        slider_ax = self.fig.add_axes([0.2, 0.01, 0.6, 0.03])
        # Create a slider to move between channels, and attach it to the slider axis we just created.
        self.slider = widgets.Slider(slider_ax, 'Channel', 0, self.avgcube.shape[0] - 1, valinit=0, valstep=1)
        # Create a button axis(or button picture) at the top right of the slice figure
        button_ax = self.fig.add_axes([0.05, 0.95, 0.06, 0.03])

        # Create a button to save the spectrum, and attach it to the button axis we just created. 
        # button_ax is the location of the button, 'Save' is the text on the button, and hovercolor is the color of the button when you hover over it.
        self.button = widgets.Button(button_ax, 'Save', hovercolor = 'gray')
        self.channel = 1
        # Attach an update function to the slider. This function will be called whenever the slider value changes.
        # What this means is that whenever you move the slider, it will call the update_channel function and update the slice.
        self.slider.on_changed(self.update_channel)
        self.ix = 0
        self.iy = 0
        self.filename= filename
    
    def update_channel(self, val):
        """ 
        Updates the slice display as you move the slider to different channels. Does not return anything, but updates the SpectrumViewer object.
        This is the function that is called when you move the slider. It takes the slider value as input, which is the channel number.

        Returns:
            None
        Parameters:
            val (int): The slider value (channel number)
        """
        
        self.channel = int(self.slider.val)
        self.display_slice(self.channel)
    
    def ticks_to_radec(self):
        """
        Convert the tick labels on the slice plot from pixel coordinates to RA and DEC.This doesn't require any input
        parameters since it uses information stored in the FITS file headers, which is already stored in the SpectrumViewer
        object in the __init__ method. Does not return anything, but updates the SpectrumViewer object.
        
        Returns:
            None
        Parameters:
            None
        """

        # Extract necessary header information for RA and DEC conversion
        # CRVAL1, CRVAL2: Reference values for RA and DEC in degrees
        # CDELT1, CDELT2: Increments per pixel for RA and DEC in degrees
        # CRPIX1, CRPIX2: Reference pixels for RA and DEC
        crval1, crval2 = self.header['CRVAL1'], self.header['CRVAL2']
        cdelt1, cdelt2 = self.header['CDELT1'], self.header['CDELT2']
        crpix1, crpix2 = self.header['CRPIX1'], self.header['CRPIX2']

        # Define the pixel coordinates where tick marks will be placed
        # For both x (RA) and y (DEC), 6 tick marks are evenly placed. Number of ticks is arbitrary and can be changed.
        x_ticks = np.linspace(0, self.avgcube.shape[2] - 1, 6)
        y_ticks = np.linspace(0, self.avgcube.shape[1] - 1, 6)

        # Convert pixel coordinates to RA and DEC using the formulas:
        # RA  = CRVAL1 + (pixel_x - CRPIX1) * CDELT1
        # DEC = CRVAL2 + (pixel_y - CRPIX2) * CDELT2
        ra_ticks = crval1 + (x_ticks - crpix1) * cdelt1
        dec_ticks = crval2 + (y_ticks - crpix2) * cdelt2
        
        # Update the tick labels on the plot to show RA and DEC values
        # Instead of showing pixel coordinates, we show the corresponding RA and DEC.
        self.slice_ax.set_xticks(x_ticks)
        self.slice_ax.set_xticklabels(ra_ticks.round(3))
        self.slice_ax.set_yticks(y_ticks)
        self.slice_ax.set_yticklabels(dec_ticks.round(3))

        # Label the x and y axes as RA and DEC respectively
        self.slice_ax.set_xlabel("RA (degrees)")
        self.slice_ax.set_ylabel("DEC (degrees)")

    def pixel_to_radec_click(self, ix, iy):
        """
        This method takes the click coordinates, taken from the capture_click_coordinates method, and converts them to RA and DEC.
        The click must be on the slice plot and is between 0 and 144, and I am doing the conversion so when someone uses this code
        and finds a source and saves the spectrum, it also saves the RA and DEC values so you can look up the coordinates elsewhere and find 
        the source in the sky.

        Returns:
            None
        Parameters:
            ix (int): The x coordinate of the click in pixel coordinates.
            iy (int): The y coordinate of the click in pixel coordinates.
        """
        # Extract necessary header information for RA and DEC conversion
        crval1, crval2 = self.header['CRVAL1'], self.header['CRVAL2']
        cdelt1, cdelt2 = self.header['CDELT1'], self.header['CDELT2']
        crpix1, crpix2 = self.header['CRPIX1'], self.header['CRPIX2']
        # Convert pixel coordinates to RA and DEC
        ra = crval1 + (ix - crpix1) * cdelt1
        dec = crval2 + (iy - crpix2) * cdelt2
        print(f"RA: {ra}, DEC: {dec}")
        return ra,dec
    
    def sum_spectrum(self):
        """
        Sum the flux values across a square aperture around the clicked point. 
        Does not require any input parameters since it uses information stored in the SpectrumViewer object in the __init__ method.
        Does not return anything, but stores the summed spectrum in the SpectrumViewer object. Does not work currently.

        Parameters:
            None

        """
        # Define an aperture square around the clicked point.
        d = 0 # Aperture size (half-size of the square)
        
        # Calculate the bounds for the square.
        x_start, x_end = max(self.ix - d, 0), min(self.ix + d, 143)
        y_start, y_end = max(self.iy - d, 0), min(self.iy + d, 143)
        print(self.avgcube[:,self.iy,self.ix])
        print('x_start',x_start,'x_end',x_end,'y_start',y_start,'y_end',y_end)
        print('ix - d = ',self.ix - d)
        # Sum over the square region across all channels.
        #summed_spectrum = np.sum(self.avgcube[:,x_start:x_end, y_start:y_end] , axis=(1, 2))
        summed_spectrum = np.sum(self.avgcube[:,x_start:x_end, y_start:y_end] , axis=(2, 1))
        print(np.shape(self.data_cube))
        self.summed_spectrum = summed_spectrum if summed_spectrum.mean() != 0 else None
        if self.summed_spectrum is not None:
            print("Shape of summed_spectrum:", self.summed_spectrum.shape,"Shape of avg_cube:", self.avgcube.shape)
   
    def display_slice(self, channel):
        """
        Display a 2D slice of the data cube for a given channel.
        
        Parameters:
            channel (int): The channel number to display.
        """
        #Clear the previous slice to prepare for the new one
        self.slice_ax.clear()

        # Display the 2D slice of the data cube for the given channel
        self.slice_ax.imshow(self.avgcube[channel, :, :], origin='lower')
        self.slice_ax.invert_xaxis()
        self.slice_ax.set_title(f"Channel: {channel}")
        self.ticks_to_radec()
        self.sum_spectrum()
        # Refresh the figure to reflect these changes
        self.fig.canvas.draw()
    

    def onclick(self, event):
        if event.inaxes == self.slice_ax:    
            """
            Handle the click event, generate the spectrum for the clicked point, and display it.
            
            Parameters:
                event (Event): Matplotlib click event.
            """

            self.capture_click_coordinates(event)
            spectrum_to_plot = self.generate_spectrum()
            self.update_spectrum_plot(spectrum_to_plot)


    def capture_click_coordinates(self, event):
        """
        Store the coordinates of the clicked point, in coordinate space. Meaning we store x and y values ranging between 0 and 144
        since that is the size of our coordinate arrays in the datacube.
        """
        self.ix, self.iy = int(event.xdata), int(event.ydata)
        print(f'User clicked at coordinates x={self.ix}, y={self.iy}')
        self.pixel_to_radec_click(self.ix, self.iy)


    def generate_frequency_array(self):
        """
        Generate a frequency array corresponding to the channels.
        
        Returns:
            numpy.ndarray: An array containing the frequency values for each channel.
            
        Formula: Frequency = CRVAL3 + (Channel - CRPIX3) * CDELT3
        """
        CRVAL3 = self.header['CRVAL3']
        CDELT3 = self.header['CDELT3']
        CRPIX3 = self.header['CRPIX3']
        
        num_channels = self.avgcube.shape[0]
        channels = np.arange(num_channels)
        
        frequencies = CRVAL3 + (channels - CRPIX3) * CDELT3
        return frequencies

    
    def generate_spectrum(self):
        """
        
        Call the summed spectrum we generate. If the summed spectrum is None, then we just plot the spectrum at the clicked point.
        If the summed spectrum is not None, then we plot the summed spectrum. Does not require any input parameters since it uses
        information stored in the SpectrumViewer object (self). 

        Returns:
            numpy.ndarray: The spectrum to plot.

        Parameters:
            None
        """
        self.sum_spectrum()
        spectrum_to_plot = self.summed_spectrum if self.summed_spectrum is not None else self.avgcube[:, self.iy, self.ix]
        self.spectrum_to_plot = spectrum_to_plot
        print(type(spectrum_to_plot))
        return spectrum_to_plot
    
    def update_spectrum_plot(self, spectrum_to_plot):
        """
        Update the spectrum plot with the new spectrum. Does not return anything, but updates the SpectrumViewer object.

        Returns:
            None
        Parameters:
            spectrum_to_plot (numpy.ndarray): The spectrum to plot.
        """
        self.spectrum_line.set_ydata(spectrum_to_plot)
        self.spectrum_line.set_xdata(range(len(spectrum_to_plot)))
        self.spectrum_line.set_color('blue')
        self.spectrum_ax.relim()
        self.spectrum_ax.autoscale_view(True, True, True)
        self.fig.canvas.draw()

    def save(self, event):
        """
        
        Save the spectrum to a FITS file. Does not return anything, but saves the spectrum to a FITS file. Takes the click event as input.
        Click event is just the click on the figure, so we can use that to check if the user clicked on the save button or not.

        Returns:
            None
        Parameters:
            event (Event): Matplotlib click event.
        """
        # Check if the user clicked on the save button. When the user clicks, that generated a click event, which just means it stores the 
        # x and y coordinates of the click. When there is click event, this if statement goes 'Was the click within the coordinates of the save button?'
        # and if it was, then it saves the spectrum.
        if event.inaxes == self.button.ax:
            print("Saving")
            
            # Get RA and DEC
            ra, dec = self.pixel_to_radec_click(self.ix, self.iy)
            
            # Generate frequency array
            frequencies = self.generate_frequency_array()
            print(f'This is the problem values', frequencies[frequencies == 0])
            velocities = ((self.header['RESTFREQ']/(frequencies*1e6) - 1)* c.to('km/s'))
            print('Velocity array is',type(velocities[0]))
            # Create a new FITS file
            tab = Table([frequencies,velocities,self.spectrum_to_plot], names = ('FREQUENCY','VELOCITY','FLUX'))
            tab.write('test_spectrum.fits', format = 'fits', overwrite = True)
            plt.savefig(f'spectrum{self.channel}_{self.ix}_{self.iy}.png',dpi = 300)
            print("Saved")
            
    def main(self, channel):
        """
        Run the SpectrumViewer for a given channel.
        
        Parameters:
            channel (int): The channel number to display initially.
        """
        # Display the initial slice
        self.display_slice(channel)
        
        '''The method mpl_connect is used to connect event-handling functions (also known as callbacks) to Matplotlib figures. 
        This is the mechanism that allows you to make your Matplotlib plots interactive. 
        The mpl_connect method belongs to the Canvas object, which is essentially the drawing area of the Matplotlib figure.'''
    
        self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        self.fig.canvas.mpl_connect('button_press_event', self.save)
        plt.show()

# Usage example
viewer = SpectrumViewer('1220+09a_spectral.fits')
viewer.main(10)
