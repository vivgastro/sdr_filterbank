from rtlsdr import RtlSdr
import numpy as np
import matplotlib.pyplot as plt
import time, sys
import argparse
from sigproc import SigprocFile as SF
import jdutil
import datetime

def get_period():
    F0 = 2.1974245233412279415
    F1 = -9.6961455414359592752e-14
    EPOCH = 57600
    EPOCH_NOW = jdutil.jd_to_mjd(jdutil.datetime_to_jd(datetime.datetime.now()))
    EPOCH_NOW = 60211
    delta_days = EPOCH_NOW - EPOCH
    delta_seconds = delta_days * 24 * 60 * 60
    delta_F = delta_seconds * F1
    F0_now = F0 + delta_F
    P0_now = 1 / F0_now
    return P0_now

def nearest_pow_2(x, max_pow = 30):
    powers_of_2 = 2**np.arange(max_pow)
    nearest = np.argmin(np.abs(x - powers_of_2))
    return powers_of_2[nearest]
    

def make_sigproc_header(cfreq, nchan, tsamp, bw, src_ra_deg=0.0, src_dec_deg=0.0, src_name = "Fake", nbits=64):
    '''
    cfreq: float, MHz
    nchan: int, count
    tsamp: float, seconds
    bw: float, MHz
    src_ra_deg: float, degrees
    src_dec_deg: float, degrees
    src_name: str, name
    nbits: int, nbits
    '''
    time_now = datetime.datetime.now()
    ch_bw = bw / nchan
    fch1 = cfreq - bw + ch_bw/2
    header = {
            'nbits':nbits,
            'nchans':nchan,
            'nifs': 1,
            'tstart': jdutil.jd_to_mjd(jdutil.datetime_to_jd(time_now)),
            'tsamp': tsamp,
            'fch1': fch1,
            'foff':ch_bw,
            'src_ra_deg': src_ra_deg,
            'src_dec_deg': sec_dec_deg,
            'source': src_name,
            }
    
    return header            


def run_filterbank(signal, tsamp, nchan):
    f, t, spec = scipy.signal.spectogram(signal, 
                                         fs = 1 / tsamp,
                                         window = scipy.signal.window.tukey(nchan, 0.25), 
                                         nperseg = nchan,
                                         noverlap = 0,
                                         nfft = nchan,
                                         scaling='spectrum',
                                         mode = 'magnitude')
    return f, t, spec


def main():

    Fs = int(args.bw * 1e6)         #Sampling rate
    Fcen = args.cfreq               #MHz
    Gain = args.gain
    nchan = args.nchan

    desired_tres = args.tres * 1e-3               #sec
    recording_len = args.dur                      #seconds
    capture_samps = 65536*8                       #largest size supported by SDR
    
    outfile = args.outname
    
    sdr = RtlSdr()
    sdr.central_freq = int(Fcen*1e6)
    sdr.sample_rate = Fs
    if Gain is not None:
        sdr.gain = Gain

    actual_tsamp = 1./sdr.sample_rate   #sec
    print("Actual Sampling rate is ", sdr.sample_rate)
    print("Actual BW is ", sdr.sample_rate * 1e-6, " MHz")
    print("Actual raw tsamp is ", actual_tsamp)


    #period = get_period()
    #capture_samps = int(period / actual_tsamp)
    print(f"Capture_samps is {capture_samps}")
    print(f"Capture_len is {capture_samps * actual_tsamp}")

    tx = nearest_pow_2(int(desired_tres / actual_tsamp / nchan ))

    print(f"tx is {tx}")
    print(f"capture_samps / tx = {capture_samps / tx}")
    print(f"final tres is {actual_tsamp * tx} s")

    output_data_rate = capture_samps / tx * nchan #per_sec
    print(f"Output date rate is {output_data_rate * 1e-3}KB/s")

    header = make_sigproc_header(Fcen, nchan, actual_tsamp, sdr.sample_rate * 1e-6)

    of = SigprocFile(args.outname, 'wb', header)
    #of = open(outfile, 'wb')

    global counter
    counter = 0
    Ncaptures = int(recording_len / (capture_samps * actual_tsamp)) - 1

    def filterbank_tscrunch_and_save(samps, contxt):
        global counter
        print(f"Count {counter}/{Ncaptures} - Got {len(samps)} samps")
        f, t, spec = run_filterbank(samps, actual_tsamp, nchan)
        fil = spec.reshape(-1, tx).sum(axis=-1)
        fil.tofile(of.fin)
        counter += 1
        if counter > Ncaptures:
            print("Closing output file")
            of.close()
            print("Closing async read")
            sdr.cancel_read_async()
            print("Exiting the script")
            sys.exit()

    sdr.read_samples_async(callback=filterbank_tscrunch_and_save, num_samples = capture_samps)
    sdr.close()


if __name__ == '__main__':
    a = argparse.ArgumentParser()
    a.add_argument("-cfreq", type=float, help="Center freq in MHz (def=820)", default=820)
    a.add_argument("-bw", type=float, help="BW in MHz (def=3.2)", default=3.2)
    a.add_argument("-dur", type=float, help="Observation duration in seconds (def=10)", default=10)
    a.add_argument("-tres", type=float, help="Desired time resoltion in ms (def=1)", default=1)
    a.add_argument("-nchan", type=int, help="Desired nchan (def=64)", default=64)
    a.add_argument("-gain", type=float, help="Desired Gain (def=auto)", default=None)
    a.add_argument("-outname", type=str, help="Name of the output file", required = True)
    args = a.parse_args()
    main()


