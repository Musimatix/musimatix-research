library(tuneR)
library(seewave)
epsilon = 0.0000001
setwd ("/Users/Kondor/Desktop/Music/FromC#/korpus")
x1 = readWave("Tolya-Rodina-pripev.wav")
x1 = updateWave(x1)

#размер окна 10мс
#шаг спектра 5мс

signal = x1@left / max(abs(x1@left))
signal = signal[(5.75*44100):(5.85*44100)]
windowlen = 0.02
steplen = windowlen/2

x = fft(signal)
a = matrix (0,length(x)/2,2)
bin = 44100/length(x)
step = 0
for (i in 1:(length(x)/2))
{
    step = step + bin
    if (step < 2000)
    {
        a[i,2] = abs(x[i])
        a[i,1] = step
    }
}
plot (a,type = 'l')



signal = x1@left / max(abs(x1@left))
signal = signal[(5.85*44100):(5.95*44100)]
windowlen = 0.02
steplen = windowlen/2

x = fft(signal)
a = matrix (0,length(x)/2,2)
bin = 44100/length(x)
step = 0
for (i in 1:(length(x)/2))
{
    step = step + bin
    if (step < 2000)
    {
        a[i,2] = abs(x[i])
        a[i,1] = step
    }
}
points (a,type = 'l',col = 'green')

signal = x1@left / max(abs(x1@left))
signal = signal[(5.95*44100):(6.05*44100)]
windowlen = 0.02
steplen = windowlen/2

x = fft(signal)
a = matrix (0,length(x)/2,2)
bin = 44100/length(x)
step = 0
for (i in 1:(length(x)/2))
{
    step = step + bin
    if (step < 2000)
    {
        a[i,2] = abs(x[i])
        a[i,1] = step
    }
}
points (a,type = 'l',col = 'red')



#------------ расширенный анализ

# выделить энергетически приятные участки, дальше работать только на них
# идея - найти ближайший пик к 500
# в случае падения амплитуды - слог




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

windowlen = 0.02
steplen = windowlen

signal = x1@left[1:(length(x1)/3)]

E = Energy(signal,windowlen*x1@samp.rate,steplen*x1@samp.rate)*100
a = matrix(0,length(signal))
k = windowlen*x1@samp.rate
for ( i in 1:length(E))
{
    a[k*i] = E[i]
    
}

plot(E,type = 'l')
#lines(a*5000, col = 'red' )
for ( i in 1:length(a))
{
    if (a[i]<0.1)
        signal[i] = 0
}
plot(signal,type = 'h')
