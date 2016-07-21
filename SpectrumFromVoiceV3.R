library(tuneR)
library(seewave)
epsilon = 0.0000001
setwd ("/Users/Kondor/Desktop/Music/FromC#")
x1 = readWave("L3.wav")
x1 = updateWave(x1)

#размер окна 10мс
#шаг спектра 5мс

signal = x1@left / max(abs(x1@left))

windowlen = 0.02
steplen = windowlen/2

delta = x1@samp.rate/length(window)

all = function (signal,wlen,step)
{
    len = length(signal)
    pos = 1
    #Число фреймов
    NumFrames = floor((len-wlen)/step) + 1
    entropy = matrix (0,NumFrames,1)
    for (j in 1: NumFrames)
    {
        window = signal[pos:(pos+wlen-1)]
        pos = pos+1
        #посмотрим на энтропию
        #вероятности
        p = matrix (0,length(window),1)
        index = 0
        for (i in 1: length(window))
        {
            index = floor(abs(window[i]*x1@samp.rate))+1
            p[index] = p[index]+1
        }
        #normilize
        p = p/sum(p)
        #entropy
        for ( i in 1:length(window))
        {
            if (p[i] > epsilon)
            {
                entropy[j] = entropy[j] + p[i] * log2 (p[i])
            }
        }
    }
}

all(x1@left,windowlen*x1@samp.rate,steplen*x1@samp.rate)
