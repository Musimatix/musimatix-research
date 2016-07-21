library(tuneR)
library(seewave)
epsilon = 0.0000001
setwd ("/Users/Kondor/Desktop/Music/FromC#")
x1 = readWave("11.wav")
x1 = updateWave(x1)

#размер окна 50мс
#шаг спектра 25мс


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
Centroid = function (signal,wlen,step)
{
    #нормируем
    signal = signal/max(abs(signal))
    len = length(signal)
    pos = 1
    m = (1:(wlen/2))/(wlen/2+1)
    #Число фреймов
    NumFrames = floor((len-wlen)/step) + 1
    Centroidf = matrix(0,NumFrames,1)
    for (i in 1:NumFrames)
    {
        #применим окно Хэмминга
        window = signal[pos:(pos+wlen-1)]
        for (j in 1:wlen)
        {
            window[j] = window[j] * (0.53836 - 0.46164 * cos (2 * pi * j / (wlen - 1)))
        }
        #Теперь Фурье (внимание! надо делить итог на 2)
        p= fft(window)
        magnitude = (Re(p)*Re(p) + Im(p)*Im(p))^0.5
        magnitude = magnitude[1:(wlen/2)]
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
steplen = windowlen/2

E = Energy(x1@left,windowlen*x1@samp.rate,steplen*x1@samp.rate)*100

#Смотрим слоги
countlocmax = 0
locmaxi = 0
steptime = matrix (0,length(E),2)
for (i in 6:length(E))
    if (((E[i]>(1.5*E[i-6])) | (E[i]>(1.5*E[i-5])) | (E[i]>(1.5*E[i-4])) | (E[i]>(1.5*E[i-3])) 
         | (E[i]>(1.5*E[i-2])) | (E[i]>(1.5*E[i-1]))) && (E[i]>E[i+1]) && (E[i]>0.8)
        && ((E[i]>(E[i-1]))))#условие на максимум
    {
        if ((i - locmaxi > 4) ) #склеим близкие
        {
            locmaxi = i
            countlocmax = countlocmax + 1
            steptime[countlocmax,2] = E[i]
            steptime[countlocmax,1] = i
        } else #переопределим близкие
        {
            if (E[i] > steptime[countlocmax,2] )
                {
                    steptime[countlocmax,2] = E[i]
                    steptime[countlocmax,1] = i
                    locmaxi = i
                }
        }
    }
plot(E,type = 'l')
points(steptime)


#steptime[1:countlocmax,] = steptime[1:countlocmax,] * steplen
#plot(x1@left[x1@left>0,type ])
#steptime[1:countlocmax,]

#windowlen = 0.01
#steplen = windowlen/2
#C = Centroid(x1@left,windowlen*x1@samp.rate,steplen*x1@samp.rate)
#plot(C,type = 'l')


#10 MEL частот(12 опрных точек)
m = c(401.25, 622.50, 843.75, 1065.00, 1286.25, 1507.50, 1728.74, 1949.99, 2171.24, 2392.49, 2613.74, 2834.99)
#10 MEL частот в герцах(12 опрных точек)
h = c(300, 517.33, 781.90, 1103.97, 1496.04, 1973.32, 2554.33, 3261.62, 4122.63, 5170.76, 6446.70, 8000)

#на спектре нашего фрейма
f = floor((windowlen*x1@samp.rate+1)*h/x1@samp.rate)

#MFCC filter

H = function (m,k,f)
{
    if (k < f[m-1])
    {
        return( 0)
    }
    if ((f[m-1]<= k) && (f[m] >= k))
    {
        return( (k - f[m-1])/(f[m]-f[m-1]))
    }
    if ((f[m]<= k) && (f[m+1] >= k))
    {
        return(  (f[m+1] - k) / (f[m+1] - f[m]))
    }
    if (k > f[m+1])
    {
        return(  0)
    }
}
#Применение фильтра
#S = matrix(0,10)
#for (j in 1:10)
#    {
#        for (i in 1:length(window))
#            S[j] = S[j]+ (abs(window[i])^2) * H((j+1),i,f)
#    S[j] = log(S[j])
#    }

MEL = function (signal,wlen,step)
{
    #нормируем
    signal = signal/max(abs(signal))
    len = length(signal)
    pos = 1
    #Число фреймов
    NumFrames = floor((len-wlen)/step) + 1
    MEL = matrix(0,NumFrames,11)
    for (k in 1:NumFrames)
    {
        #применим окно Хэмминга
        window = signal[pos:(pos+wlen-1)]
        pos = pos + 1
        S = matrix(0,10)
        for (j in 1:10)
        {
            for (i in 1:length(window))
                S[j] = S[j]+ (abs(window[i])^2) * H((j+1),i,f)
            #if (S[j] > epsilon)
            #    S[j] = log(S[j])
        }
        MEL[k,(2:11)] = S
        MEL[k,1] = sum (S^2)
    }
    return(MEL)
}

#-----------------------------
bpm = 180
plot(E,type = 'l')
steptime = steptime[1:countlocmax,]
points(steptime)
timestep = steptime[1:countlocmax,] * steplen
#ударов всего около (floor(timestep[countlocmax] * bmp / 60))
delta = length(E)/(floor(timestep[countlocmax] *bpm / 60))
beat = matrix(0,(floor(timestep[countlocmax] * bpm / 60)+1),2)
beat[1,1] = 0
beatmax = 0
beatmaxt = 0
for ( i in 2:length(beat[,1]))
{
    beat[i,1] = beat[i-1,1] + delta
}
beat[,2] = 0

for (j in 0:(floor(delta)-1))
{
    a = beat[,1] + j
    b = E[floor(a[-length(a)])]
    b = b[b>(max(b)/2)]
    if (sum(b) > beatmax)
    {
        beatmaxt = j
        beatmax = sum(E[beat[,1]])
    }
}
beat[,1] = beat[,1] + beatmaxt
points(beat,col = 'red')



#-----------------------------
min = length(E)
for ( i in 1:(countlocmax-1))
{
    for (j in (i+1):countlocmax)
        if (abs(steptime[i,1] -steptime[j,1]) < min)
        {
            min = abs(steptime[i,1] - steptime[j,1])
        }
}
min
beatmin = matrix(0,(floor(length(E) / min )),2)
beatmin[1,1] = 0
for ( i in 2:(length(beatmin[,1])-1))
{
    beatmin[i,1] = beatmin[i-1,1] + min
}
beatmin[,2] = 1
beatmax = 0
beatmaxt = 0
for (j in 0:(floor(min)-1))
{
    a = beatmin[,1] + j
    b = E[floor(a[-length(a)])]
    b = b[b>(max(b)/2)]
    if (sum(b) > beatmax)
    {
        beatmaxt = j
        beatmax = sum(E[floor(a[-length(a)])])
    }
}
beatmin[,1] = beatmin[,1] + beatmaxt

points(beatmin,col = 'green',type = 'h')

#--------------------------
words = matrix(0,countlocmax)
begin = 0
step = 1
#Выделение слов
for ( i in 1:length(E))
{
    if ((begin == 0) && (E[i]>0.25)) #начало возможного слова
    {
        begin = i
    }
    if ((begin != 0) && (i-begin > 30) && (E[i]<0.25)) #end of the word
    {
        words[step] = begin
        words[step+1] = i-1
        begin = 0
        step = step + 2
    }
    if ((begin != 0) && (i-begin < 30) && (E[i]<0.25)) #too short word
        begin = 0
}
#points(words)
wordf = matrix(0,length(words),2)
wordf[,1] = words
points(wordf)
words = words[1:(step-1)]
words
#Анализ слов по отдельности
if (sum(words) != 0){
for ( i in 1:(length(words)/2))
{
    k = i*2
    d = matrix(0,(words[k]-words[k-1]-1))
    step = 1
    for (j in (words[k-1]+1):words[k])
    {
        d[step] = (E[j]-E[j-1])/2
        step = step+1
    }
    #plot(d,type = 'h')
    
}}
