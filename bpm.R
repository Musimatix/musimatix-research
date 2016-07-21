library(tuneR)
library(seewave)
epsilon = 0.0000001
setwd ("/Users/Kondor/Desktop/Music/FromC#")
x1 = readWave("USSR.wav")
x1 = updateWave(x1)


BPM = function (signal, AmpMax, wlen = 5*44100, samprate = 44100 , beg = 120, rate = 60)
{
    #нормируем
    signal = signal/max(signal)
    ta = matrix (0,wlen)
    tb = matrix (0,wlen)
    window = signal[1:wlen]
    for (j in 1:wlen)
    {
        window[j] = window[j] * (0.53836 - 0.46164 * cos (2 * pi * j / (wlen - 1)))
    }
    FFT = fft(window)
    #Test bpm
    goodchoiceb = 0
    goodchoicee = 0
    koef = rate/6
    E = matrix(0, length((beg-rate)/koef):((beg+rate)/koef))
    for (i in ((beg-rate)/koef):((beg+rate)/koef))
    {
        bpm = i*koef
        Ti = 60/bpm *samprate    
        l = matrix(0,wlen)
        j = matrix(0,wlen)
        for (k in 1:wlen)
        {
            if ((k %% Ti) == 0)
            {
                l[k] = 1
            }
        }
        E[i] = sum (fft(l) * FFT)
        e = abs(E[i])
        if(e > goodchoicee)
        {
            goodchoicee = e
            goodchoiceb = bpm
        }
    }
    return (goodchoiceb)
}

Sys.time()
bpm = BPM(x1@left[(length(x1@left)/2 - 2.5 * 44100 ):(length(x1@left)/2 + 2.5 * 44100 )], 
          wlen = 5*44100, AmpMax = 32767, 
          samprate = x1@samp.rate , beg = 120, rate = 60 )
bpm1 = BPM(x1@left[(length(x1@left)/2 - 1.5 * 44100 ):(length(x1@left)/2 + 1.5 * 44100 )], 
          wlen = 3*44100, AmpMax = 32767, 
          samprate = x1@samp.rate , beg = bpm, rate = 12 )

#bpm2 = BPM(x1@left[1:(5 * 44100)], wlen = 5*44100, AmpMax = 32767, 
#          samprate = x1@samp.rate , beg = bpm1, rate = 6 )
Sys.time()
bpm
#bpm1
#bpm2
