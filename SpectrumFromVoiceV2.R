library(tuneR)
library(seewave)
setwd ("/Users/Kondor/Desktop/Music/FromC#")
x1 = readWave("11.wav")
x1 = updateWave(x1)

#размер окна 50мс
#шаг спектра 50мс
#интервал фрейма 40мс

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

#Спектральный центроид
#signal - исходный сигнал
#wlen - длина окна
#step - шаг
#E - энегрия
Centroid = function (signal,wlen,step,E)
{
    #нормируем
    signal = signal/max(signal)
    len = length(signal)
    pos = 1
    m = (2:(wlen/2+1))/(wlen/2)
    #Число фреймов
    NumFrames = floor((len-wlen)/step) + 1
    Centroidf = matrix(0,NumFrames,1)
    for (i in 1:NumFrames)
    {
        #применим окно Хэмминга
        window = signal[pos:(pos+wlen-1)]
        for (j in 1:wlen)
        {
            window[j] = window[j] * (0.53836 - 0.46164 * cos (2 * pi * j / (len - 1)))
        }
        #Теперь Фурье (внимание! надо делить итог на 2)
        magnitude = abs(fft(window)[1:(wlen/2)])
        magnitude = magnitude/max(magnitude)
        Centroidf[i] = sum (m * magnitude)/ sum (magnitude)
        #a = meanspec(window,44100,plot = F)
        #Centroidf[i] = specprop(a,44100)$cent
        
        #Если энергия мала - это просто шум
        if (sum(window^2) < 0.01){
            Centroidf[i] = 0}
        pos = pos + step
        #print (Centroidf[i])
    }
    Centroidf
}

windowlen = 0.05
steplen = 0.02

E = Energy(x1@left,windowlen*x1@samp.rate,steplen*x1@samp.rate)*100
plot(E,type = 'l')
#C = Centroid(x1@left,windowlen*x1@samp.rate,steplen*x1@samp.rate,E)
#plot(C,type = 'l')


#Смотрим слоги
countlocmax = 0
locmaxi = 0
steptime = matrix (0,length(E),2)
for (i in 6:length(E))
    if ((E[i]>(E[i-3])) && (E[i]>E[i+1]) && (E[i]>1)) #условие на максимум
    {
        if ((i - locmaxi)>5) #склеим близкие
        {
            locmaxi = i
            countlocmax = countlocmax + 1
            steptime[countlocmax,2] = E[i]
            steptime[countlocmax,1] = i#*windowlen*length(x1@left)/(length(E))/x1@samp.rate
        }
    }
steptime[1:countlocmax,]
points(steptime[1:countlocmax,])

