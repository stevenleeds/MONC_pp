"""
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Copyright (C) 2012 Juerg Schmidli, ETHZ

=================
SPECTRA package
=================

Collection of routines for calculation of turbulence spectra.

References:
* Chapter 7 in Kaimal and Finnigan, 1994
* Section 5.2 in Mazzoni, 1996 (PhD thesis)
"""

import numpy as np
import numpy.fft as fft
 
# -----------------------------------------------------------------------------
def fft_plus(data, Fs=1.):
    """ Compute the one-sided power spectral density.

    Compute the one-sided power spectral density using FFT. The spectrum is
    normalized such that sum(p) = var(data).

    Parameters:
    -----------
    data : array_like
        Array containing the data.
    Fs : scalar, optional
        The sampling frequency.

    Returns:
    --------
    (p, freq) : ndarray, ndarray
        
    Refs:
        1702 - Using FFT to Obtain Simple Spectral Analysis Plots
        see http://www.mathworks.ch/support/tech-notes/1700/1702.html
    """

    n = len(data)
    p = fft.fft(data)
    nUniquePts = np.ceil((n+1)/2.0)
    p = p[0:nUniquePts]
    p = np.abs(p)
    p = p / np.float(n)
    p = p**2                # square to get power 

    # multiply by two (since half of points were dropped)
    if n % 2 > 0:           # we've got odd number of points fft
        p[1:len(p)] = p[1:len(p)]*2
    else:
        p[1:len(p)-1] = p[1:len(p)-1]*2     # we've got even number...
    freq = np.arange(0, nUniquePts, 1.0) * (np.float(Fs) / n);
    return (p, freq)


# -----------------------------------------------------------------------------
def window_cosine(data, percent=10):
    """ Apply Cosine-Taper to the data (to the first and last decimal of the array) """

    data = np.atleast_1d(data)
    if data.size == 1:
        window_size = (100-percent)/100.
        sfact = 1.0/window_size
        return sfact     
    N = data.size
    l = np.int(data.size*percent/100.)

    wgt = np.ones(N)
    idx = np.arange(l)
    wgt[idx] = 0.5*(1-np.cos(np.pi*idx/(N-1)*(N-1.)/l))
    idxr = N-idx-1
    wgt[idxr] = wgt[idx]
    res = wgt*data
    return res


# -----------------------------------------------------------------------------
def spectrum_basic(data, Fs=1., pad=True, detrend=None, window=window_cosine, 
        rmzf=True, scale_by_freq=False):
    """ Compute the one-sided power spectral density.

    Compute the one-sided power spectral density using FFT including the
    following optional data handling steps:
    (1) detrend the data
    (2) padding (addition of zeros)
    (3) tapering the data (tapering window) and remove non-zero mean
    with a scaling to correct for power loss due to padding and tapering. 
    
    Normalization of the return spectrum is as follows:
    (a) scale_by_freq=False:
        var(data) = sum(p)
    (b) scale_by_freq=True:
        var(data) = sum(df*p)

    Parameters:
    -----------
    data : array_like
        Array containing the data.

    Fs : scalar, optional
        The sampling frequency (samples per time/space unit). It is used to 
        calculate the Fourier frequencies, freq, in cycles per time/space
        unit. The default value is 1.

    pad : boolean
        Specifies whether the data series should be padded with zeros in order
        to obtain a block size of 2**n for the FFT.

    detrend : callable
        The function applied to each segment before the FFT, designed to
        remove the mean or linear trend.

    window : callable
        The function applied to window the data series.

    rmzf : boolean
        Specifies whether the zero-frequency should be removed from the output.

    scale_by_freq : boolean
        Specifies whether the resulting density values should be scaled by
        the scaling frequency, which gives density in units of Hz^-1.

    Returns:
    --------
    (p, freq) : ndarray, ndarray
        
    Refs and Notes:
    ---------------
    * Kaimal and Finnigan, 1994 (Chapter 7)
    * Trend removal should be performed only if trends are physically expected
        or clearly apparent in the time series. Automatic detrending is not
        recommended, except for certain variance and flux calculations (e.g. u'T')
        (KF94, page 264).
    * It is essential that means and trends be removed from the data before
        adding zeros (KF94, page 268)
    * After windowing, the finite data series may have aquired a non-zero mean
        which must be removed before applying the FFT (KF, page 267)
    """

    nt = data.size
    df = Fs/float(nt)   # width of elementary frequency band
    nfft = 2**np.int(np.ceil(np.log(data.size)/np.log(2)))

    # pad with zeros, if required
    if nfft > nt and pad:
        N_z = nfft - nt
        if N_z > nfft/3:
            print 'Warning: large padding fraction N_z > N/3!'
        n1 = N_z//2
        n2 = N_z - n1
        data2 = np.zeros(nfft, dtype=data.dtype)
        data2[n1:nfft-n2] = data
        data = data2
        sfact = nfft/np.float(nt)
    else:
        nfft = nt
        sfact = 1.0

    # tapering
    if window is not None:
        data = window(data)
        sfact_tap = window(nfft)
        sfact = sfact*sfact_tap
        # after tapering, data may have aquired a non-zero mean, remove 
        data = data - data.mean()

    # FFT
    (p, freq) = fft_plus(data, Fs=Fs)

    # correct for power loss due to padding and tapering
    p = sfact*p

    # scale by frequency
    if scale_by_freq:
        p = p/df

    if rmzf:
        return (p[1:], freq[1:])
    else:
        return (p, freq)


# -----------------------------------------------------------------------------
def spectrum(data, Fs=1., mode='def', nrec=8, smooth=False, rm_first=0, **keywords):
    """ Compute the one-sided power spectral density.

    Compute the one-sided power spectral density including optional frequency
    smoothing and spectral splicing. Resulting density values are scaled by
    the scaling frequency, such that var(data) = sum(p*df). The equality is 
    only approximate due to detrending, padding and tapering.

    The following options for computing the spectra exist:
    
    mode='def' (default):
        A single spectrum is calculated for the data series.

    mode='low':
        Compute a low-frequency spectrum by block-averaging the original data
        series prior to the spectral analysis. 

    mode='high':
        Compute a high-frequency spectrum by averaging succesive spectra
        from a set of short records obtained by subdividing the data series
        into non-overlaping blocks.
         
    mode='splice':
        Splice a low-frequency and high-frequency spectrum. Frequency smoothing
        is always used (parameter smooth is ignored).

    Parameters:
    -----------
    data : array_like
        Array containing the data.

    Fs : scalar
        The sampling frequency (samples per time/space unit). It is used to 
        calculate the Fourier frequencies, freq, in cycles per time/space
        unit. The default value is 1.

    mode : string
        Specifies type of spectra to be computed.

    nrec : integer
        Block size for block-averaging when computing low-frequency spectra.
        Number of blocks (short records) to split the data for computing high-
        frequency spectra.

    smooth : boolean
        Specifies whether frequency smoothing should be applied to the spectra.

    rm_first : integer
        Number of low-frequency bins to be removed from the output.

    Returns:
    --------
    (p, freq, df) : ndarray, ndarray, (scalar or ndarray)

    Refs and Notes:
    ---------------
    * Kaimal and Finnigan, 1994 (Chapter 7)
    * Mazzoni, 1996 (Section 5.2)
    """ 

    df = Fs/np.float(data.size)
    if mode == 'def':
        df = Fs/np.float(data.size)
        # default spectra
        (p, freq) = spectrum_basic(data, Fs=Fs, scale_by_freq=True, **keywords)

    if mode == 'splice':
        smooth = True
        data1 = data.copy()
        data2 = data.copy()
    else:
        data1 = data
        data2 = data

    if mode == 'low' or mode == 'splice':
        # calculate low-frequency spectrum: block-average the original series
        nt0 = data1.size
        nt = nt0//nrec
        if nt0 > nt*nrec:
            data1 = data1[0:nt*nrec]
        data1.shape = (nt, nrec)
        data1 = data1.mean(axis=1)
        (p, freq) = spectrum_basic(data1, Fs=Fs/np.float(nrec), scale_by_freq=True, **keywords)
        if mode == 'splice':
            (p_low, freq_low) = (p, freq)

    if mode == 'high' or mode == 'splice': 
        # calculate a high-frequency spectrum: average spectra over several short periods
        nt0 = data2.size
        nt = nt0//nrec
        if nt0 > nt*nrec:
            data2 = data2[0:nt*nrec]
        data2.shape = (nrec, nt)
        p = 0
        df = Fs/np.float(nt)
        for k in xrange(nrec):
            (ptmp, freq) = spectrum_basic(data2[k,:], Fs=Fs, scale_by_freq=True, **keywords)
            p = p+ptmp
        p = p/nrec
        if mode == 'splice':
            (p_high, freq_high) = (p, freq)

    if smooth:
        if mode == 'splice':
            (p_low, freq_low, df_low) = freq_smoothing(p_low, freq_low)
            (p_high, freq_high, df_high) = freq_smoothing(p_high, freq_high)
            # search for splice frequency
            n1 = 1
            n2 = 1  
            while True:
                if freq_low[-n1-1] < freq_high[n2]:
                    break
                else:
                    n2 = n2+1
                    if freq_low[-n1-1] < freq_high[n2]:
                        break
                    else:
                        n1 = n1+1
            #print n1, n2, freq_low[-n1-1], freq_high[n2]
            p_low = p_low[:-n1]
            freq_low = freq_low[:-n1]
            df_low = df_low[:-n1]
            p_high = p_high[n2:]
            freq_high = freq_high[n2:]
            df_high = df_high[n2:]
            p = np.concatenate((p_low, p_high))
            freq = np.concatenate((freq_low, freq_high))
            df = np.concatenate((df_low, df_high))
        else:
            (p, freq, df) = freq_smoothing(p, freq)

    # remove lowest frequencies
    if rm_first > 0:
        freq = freq[rm_first:]
        p = p[rm_first:]

    return (p, freq, df) 
    

# -----------------------------------------------------------------------------
def spectrum_comp(data, nt_block=None, norm=None, **keywords):
    """ Calculate power spectral density """

    if nt_block is None:
        nt_block = 1024

    nblocks = np.int(np.floor(data.size//nt_block))
    if nblocks == 0:
        raise ValueError('data array too small for nt_block = ' + str(nt_block))

    if nt_block == 1:
        (p, freq, df) = spectrum(data, **keywords)
        p = p/norm
    else:
        if nt_block*nblocks < data.size:
            data = data[:nt_block*nblocks]
        data.shape = (nblocks, nt_block)
        pm = 0
        for n in xrange(nblocks):
            (p, freq, df) = spectrum(data[n,:], **keywords)
            if norm is not None:
                p = p/norm[n]
            pm = pm + p
        p = pm/nblocks 

    return (p, freq, df)


# -----------------------------------------------------------------------------
def freq_smoothing(p, freq, nbnds_only=False):
    """ Implement frequency smoothing after Mazzoni 1996 (page 52),
        based on Kaimal and Gaynor (1983) """

    bno_256 = np.array((1,1,2,2,3,3,4,4,5,6,7,8,9,11,12,15,16,20,21,26,\
            27,34,35,45,46,59,60,77,78,100,101,130,131,169,170,220,221,256))
    bno_512 = np.array((1,1,2,2,3,3,4,4,5,5,6,7,8,9,10,12,13,16,17,21,22,28,29,37,\
        38,48,49,63,64,82,83,107,108,139,140,181,182,236,237,307,308,399,400,512))
    s_fact = 0.26

    if p.size == 256:
        bno = bno_256
    elif p.size == 512:
        bno = bno_512
    else:
        nfreq = freq.size
        nmax = np.ceil(np.log10(p.size))*10
        df = np.ones(nmax)
        df0 = np.array((0.5,1,1,1,1,1,2,2,3,4,5))
        df[:df0.size] = df0
        fc_high = df0.sum()
        count = df0.size
        exit_loop = False
        while True:
            df1 = np.int(np.round(s_fact*fc_high))
            if fc_high+df1-0.6 > nfreq:
                fc_high2 = nfreq+0.5
                df1 = fc_high2-fc_high
                exit_loop = True
            fc_high = fc_high + df1
            df[count] = df1
            count = count+1
            if exit_loop: break 
        df = df[:count]
        # check that last frequency band is not too small
        if df[-1] < 0.5*df[-2]:
            df[-2] = df[-2]+df[-1]
            df = df[:-1]
        fc_bnd = df.cumsum()
        df = df[1:]
        bno = np.zeros((df.size,2))
        bno[:,0] = fc_bnd[:-1]
        bno[:,1] = fc_bnd[1:]
                
    bno.shape = (-1,2)
    bno = bno-1
    nbnd = bno.shape[0]
    if nbnds_only:
        return nbnd
    pc = np.zeros((nbnd,))
    fc = np.zeros((nbnd,))
    df = np.zeros((nbnd,))
    df0 = freq[0]
    for k in xrange(nbnd):
        pc[k] = p[bno[k,0]:bno[k,1]+1].mean()
        fc[k] = freq[bno[k,0]:bno[k,1]+1].mean()
        df[k] = freq[bno[k,1]] - freq[bno[k,0]] + df0

    return (pc, fc, df)


# -----------------------------------------------------------------------------
def nbnds_freq_smoothing(nfft):
    """ Return number of frequency bands after frequency smoothing.
    """

    tmp = np.arange(nfft)+1
    nbnds = freq_smoothing(tmp, tmp, nbnds_only=True)
    return nbnds


# -----------------------------------------------------------------------------
def spectrum_peri(data, Fs=1., rmzf=True, scale_by_freq=True, smooth=True, pad=False):
    """ Compute the one-sided power spectral density.

    Compute the one-sided power spectral density using FFT. 
    
    Normalization of the return spectrum is as follows:
    (a) scale_by_freq=False:
        var(data) = sum(p)
    (b) scale_by_freq=True:
        var(data) = sum(df*p)

    Parameters:
    -----------
    data : array_like
        Array containing the data.

    Fs : scalar, optional
        The sampling frequency (samples per time/space unit). It is used to 
        calculate the Fourier frequencies, freq, in cycles per time/space
        unit. The default value is 1.

    rmzf : boolean
        Specifies whether the zero-frequency should be removed from the output.

    scale_by_freq : boolean
        Specifies whether the resulting density values should be scaled by
        the scaling frequency, which gives density in units of Hz^-1.

    smooth : boolean
        Specifies whether frequency smoothing should be applied to the spectra.

    pad : boolean
        Specifies whether the data series should be padded with zeros in order
        to obtain a block size of 2**n for the FFT.

    Returns:
    --------
    (p, freq) : ndarray, ndarray
        
    Refs and Notes:
    ---------------
    For more details see spectrum_basic
    """

    data = np.atleast_2d(data)
    (nk, nt) = data.shape
    df = Fs/float(nt)   # width of elementary frequency band

    # pad with zeros, if required
    nfft = 2**np.int(np.ceil(np.log(nt)/np.log(2)))
    if nfft > nt and pad:
        N_z = nfft - nt
        if N_z > nfft/3:
            print 'Warning: large padding fraction N_z > N/3!'
        n1 = N_z//2
        n2 = N_z - n1
        data2 = np.zeros((nk, nfft), dtype=data.dtype)
        data2[:, n1:nfft-n2] = data
        data = data2
        sfact = nfft/np.float(nt)
    else:
        nfft = nt
        sfact = 1.0

    # tapering
    if pad:
        for k in xrange(nk):
            data[k,:] = window_cosine(data[k,:])
        sfact_tap = window_cosine(nfft)
        sfact = sfact*sfact_tap
        # after tapering, data may have aquired a non-zero mean, remove 
        data = data - data.mean(axis=1)[:,np.newaxis]

    # calc FFT
    p = 0.
    for k in xrange(nk):
        # FFT
        (ptmp, freq) = fft_plus(data[k,:], Fs=Fs)
        p = p + ptmp
    p = p/nk

    # correct for power loss due to padding and tapering
    p = sfact*p

    # scale by frequency
    if scale_by_freq:
        p = p/df

    # remove zero-frequency
    if rmzf:
        p = p[1:]
        freq = freq[1:]

    # frequency smoothing
    if smooth:
        (p, freq, df) = freq_smoothing(p, freq)

    return (p, freq)
