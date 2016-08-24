library(tuneR)
library(seewave)
epsilon = 0.0000001
setwd ("/Users/Kondor/Desktop/Music/FromC#/korpus")
x1 = readWave("Tolya-Rodina-pripev.wav")
x1 = updateWave(x1)

#Посчитаем энергию сигнала
#signal - исходный сигнал
#wlen - длина окна
#step - шаг
Energy = function (signal, wlen, step)
{
    #нормируем
    signal = signal/max(signal)
    len = length(signal)
    pos = 1
    #Число фреймов
    NumFrames = floor((len-wlen)/step) + 1
    Energyf = matrix(0,NumFrames,1)
    for (i in 1:NumFrames)
    {
        window = signal[pos:(pos+wlen-1)]
        Energyf[i] = 1 / wlen * sum (abs(window)^2) #????????? может, в квадрате
        pos = pos + step
    }
    Energyf
}

windowlen = 0.05 #(mc)
steplen = windowlen # lenght of step

E = Energy(x1@left,windowlen*x1@samp.rate,steplen*x1@samp.rate)*100 #*100
Ewindowlen = windowlen*x1@samp.rate #lenght of window, which we use in Energy detection

#---------------------------------------
# now we want to look at every word

words = matrix(0,countlocmax)
begin = 0
step = 1
#Выделение слов
for ( i in 1:length(E))
{
    if ((begin == 0) && (E[i]>0.1)) #начало возможного слова
    {
        if (step != 1)
        {
            if (i - words[step-1] > 5) #previous end too close to new begin
            {
                begin = i
            } else
            {
                words[step-1] = i
            }
        } else
        {
            begin = i
        }
    }
    if ((begin != 0) && (i-begin > 10) && (E[i]<0.1)) #end of the word
    {
        words[step] = begin
        words[step+1] = i-1
        begin = 0
        step = step + 2
    }
    #if ((begin != 0) && (i-begin < 10) && (E[i]<0.1)) #too short word
    #    begin = 0
}
words = words[words != 0]
wordsh = matrix(-20000,length(words),2)
wordsh[,1] = words*windowlen*x1@samp.rate
#points(wordf,type = 'h')
words

#generate points of local min and max of signal (fhi)
minmax = function (fhi)
{
    minimax = matrix(1,(length(fh)/10))
    steps = 1
    #распознавание последнего пика так как после может не быть спада спада
    fh = rep (0,(length(fh)+15))
    for (i in 1:(length(fhi)))
        fh[i] = fhi[i]
    # Условие на макс
    for (i in 6:(length(fhi)-5))
    {
        if (fh[i]>fh[i-1] &&fh[i]>fh[i-2] &&fh[i]>fh[i-3] &&fh[i]>fh[i-4] &&fh[i]>fh[i-5] &&
            fh[i]>fh[i+1] &&fh[i]>fh[i+2] &&fh[i]>fh[i+3] &&fh[i]>fh[i+4] &&fh[i]>fh[i+5] )
        {
            minimax[steps] = i
            steps = steps + 1
        }
    }
    # Условие на мин
    for (i in 6:(length(fhi)-5))
    {
        if (fh[i]<= fh[i-1] &&fh[i] <= fh[i-2] &&fh[i] <= fh[i-3] &&fh[i] <= fh[i-4] &&fh[i] <= fh[i-5] &&
            fh[i] <= fh[i+1] &&fh[i] <= fh[i+2] &&fh[i] <= fh[i+3] &&fh[i] <= fh[i+4] &&fh[i] <= fh[i+5] 
            && (fh[i+1]>epsilon) #последнее условие из-за большого количества нулей
        )
        {
            minimax[steps] = i
            steps = steps + 1
        }
    }
    return (minimax)
}

