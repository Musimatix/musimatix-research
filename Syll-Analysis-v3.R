#Изменение - другой алгоритм выявления слогов
library(tuneR)
library(seewave)
library(TTR)
library(wavelets)
setwd ("/Users/Kondor/Git/musimatix-research")
epsilon = 0.0000001
setwd ("/Users/Kondor/Desktop/Music/FromC#/korpus")
x1 = readWave("Tolya-Rodina-kuplet.wav")
x1 = updateWave(x1)
signal = x1@left
samprate = x1@samp.rate


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
    Energyf = Energyf/sum(Energyf)
    return (Energyf)
}


# b - параметр SMA
meanEnergy = function (signal, wlen, step, b)
{
    EN = Energy(signal, wlen, step)
    meanEnergyf = rep(0,length(EN))
    meanEnergyf[(b/2+1):(length(EN)-b/2)] = SMA(EN,b)[(b+1):length(EN)]
    meanEnergyf[(length(EN)-b/2+1):length(EN)] = rep(0,b/2)
    meanEnergyf[1:(b/2)] = rep(0,b/2)
    #for (i in (b/2+1):(length(EN)-b/2))
    #{
    #    meanEnergyf[i] = mean(EN[(i-b/2):(i+b/2)])
    #}
    return (meanEnergyf)
}


#теперь найдем слова

windowlen = 0.05 #(mc)
steplen = windowlen # lenght of step

E = Energy(signal,windowlen*samprate,steplen*samprate)*100 #*100
Ewindowlen = windowlen*samprate #lenght of window, which we use in Energy detection

#---------------------------------------
# now we want to look at every word

words = matrix(0,500)
begin = 0
step = 1
#Выделение слов
for ( i in 1:length(E))
{
    if ((begin == 0) && (E[i]>0.01)) #начало возможного слова
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
    if ((begin != 0) && (i-begin > 10) && (E[i]<0.01)) #end of the word
    {
        words[step] = begin
        words[step+1] = i-1
        begin = 0
        step = step + 2
    }
}
words = words[words != 0]
wordsh = matrix(-20000,length(words),2)
wordsh[,1] = words*windowlen*samprate


#Рассматриваем только участки с большой энергией, что бы не было шумов

ME = mean(E)/1.2 #!!!!!!!!!!!!!!!!!!!  подобрать коэфф!!
for ( i in 1:(length(signal)-Ewindowlen))
{
    #print (E[floor(i/Ewindowlen)+1])
    if (E[floor(i/Ewindowlen)+1]< 0.01) 
    {
        signal[i] = 0
    }
}

#Теперь с нахлестом!
windowlen = 0.05 #(mc)
steplen = windowlen/2 # lenght of step

#Бегаем по каждому слову и ищем слоги
pause = matrix(30000,1000,2) #подумать!!!!!!!!!!!!!!!!!!!!!!!!
t = 1 
E = Energy(signal,windowlen*samprate,steplen*samprate)
ME = meanEnergy(signal,windowlen*samprate,steplen*samprate,50)
flag = 1 #сейчас ищем начало
for (j in 2:(length(ME)-1))
{
    #начала слогов
    if (E[j]>ME[j] && E[j-1]<=ME[j-1] && flag == 1)
    {
        pause[t,1] = (j+1)*steplen*samprate
        t  = t +1 
        flag = 0
    }
    #Концы слогов
    if (E[j]>ME[j] && E[j+1]<ME[j+1] && flag == 0)
    {
        pause[t,1] = (j+1)*steplen*samprate
        t  = t +1 
        flag = 1
    }
}

#Проверка
plot(E[1:1000],ylim = c(-0.002,0.002) , type = 'h')
lines(ME,type = 'l',col = 'red')
pause1 = pause
pause1[,1] = pause[,1]/steplen/samprate - 1
pause1[,2] = -pause[,2]
lines(pause1,type = 'h',col = 'red')

plot(x1@left[1:500000],type = 'l')
points(wordsh,col = 'green',type = 'h')
points(pause,type = 'h',col = 'red')


plot(x1@left[(words[12]*samprate*windowlen):(words[13]*samprate*windowlen)],type = 'l')
points(wordsh,col = 'green',type = 'h')
points(pause,type = 'h',col = 'red')

#----------------------------------

#на основе длин и громкости внутри каждого слова
#--------------

allslen = matrix(0,nrow(pause))
stepslen = 2
allslen[1] = 1
allsvol = matrix(0,nrow(pause))
stepsvol = 2
allsvol[1] = 1
allsmvol = matrix(0,nrow(pause))
stepsmvol = 2
allsmvol[1] = 1

for (j in 1:(length(words)/2))
#for (j in 1:1)
{
    #извлекаем нужные слоги
    #начало первого слога
    sbeg = nrow(pause)
    for (i in nrow(pause):1)
    {
        if (pause[i,1]>=(words[2*j-1]*windowlen*samprate))
            sbeg = i
    }
    #Конец последнего
    send = nrow(pause)
    for (i in 1:nrow(pause))
    {
        if (pause[i,1]<=(words[2*j]*windowlen*samprate) && pause[i,1] != pause[nrow(pause),1])
            send = i
    }
    #длины слогов
    slen = rep(0,(send-sbeg))
    #мксимальная громкость слогов
    smvol = rep(0,send-sbeg)
    #"интеграл" слогов
    svol = rep(0,send-sbeg)
    #Бегаем по слогам
    for (i in sbeg:(send-1))
    {
        #берем только ненулевые значения сигнала from pause[i,1] to pause[i+1,1]
        for (l in pause[i,1]:pause[i+1,1] )
        {
            if (x1@left[l]>epsilon)
            {
                slen[i-sbeg+1] = slen[i-sbeg+1] + 1
            }
        }
        #берем среднее по 10 самым большим значениям кроме первого
        smvol[i-sbeg+1] = mean (pmax(abs(x1@left[pause[i,1]:pause[i+1,1]]))[2:10])
        #Интегральная сумма
        svol[i-sbeg+1] = sum(abs(x1@left[pause[i,1]:pause[i+1,1]]))
    }
    #нормировка
    slen = slen/sum(slen)
    smvol = smvol/sum(smvol)
    svol = svol/sum(svol)
    
    allslen[stepslen:(stepslen+length(slen)-1)] = slen
    allslen[stepslen+length(slen)] = j+1
    stepslen = stepslen+length(slen)+1
    
    allsvol[stepsvol:(stepsmvol+length(svol)-1)] = svol
    allsvol[stepsvol+length(svol)] = j+1
    stepsvol = stepsvol+length(svol)+1
    
    allsmvol[stepsmvol:(stepsmvol+length(smvol)-1)] = smvol
    allsmvol[stepsmvol+length(smvol)] = j+1
    stepsmvol = stepsmvol+length(smvol)+1
    
}

#результат - таблица, по стобцам параметры, по строчкам слоги, перед каждым словом стоит его номер
allpar = data.frame (allslen,allsvol,allsmvol)

end = 1
for (i in (nrow(allpar)-1):1)
{
    if ((allpar[i,1] == 0) && (allpar[i,2] == 0) && (allpar[i,3] == 0) && (allpar[i+1,1] == 0) && (allpar[i+1,2] == 0) && (allpar[i+1,3] == 0))
        end = i
}
allpar = allpar [1:end,]
head(allpar,150)


#хуйня какая-то

windowlen = 0.02
steplen = windowlen/2
x = x1@left[1:150000]

E = Energy(x,windowlen*samprate,steplen*samprate)

#for (i in 1:length(signal))
#избавимся от нулей
k = 0
for(i in 1:length(E))
{
    if(E[i] > 0)
        k = k+1
    #if (E[i] == 0)
    #    E[i] = epsilon
}
E1 = rep(0,k)
l = 1
for(i in 1:length(E))
{
    if (E[i] > 0)
    {
        E1[l] = E[i]
        l = l+1
    }
}
E=E1*1000
plot(E,type ='l')

#симметризуем
simE = rep(0,(2*length(E)))
simE[1:length(E)] = E
simE[length(simE):(length(E)+1)] = E
plot(simE,type ='l')

#Обратная функция и возведение в степень
revE = simE^(0.01)
revE = 1/revE
plot(revE,type = 'l')

#Коэффициенты
c = fft(revE,inverse = T)/length(revE)
c = c[Im(c)>0]
#c[1:(length(c)/2-1)] = c[2:(length(c)/2)]
#c[(length(c)/2):(length(c))] = 0
plot(c,type = 'l')

plot(unwrap(c),type = 'l')
#phaseplot(abs(c),x1@samp.rate)
#cepstro(c,x1@samp.rate,wl=x1@samp.rate*0.02)


psi = -atan2(Im(c),Re(c))
ts = seq(0,pi,pi/99)

psi = c
#производная
psi1 = rep(0,length(psi))
for (i in 2:length(psi))
{
    psi1[i] = (psi[i]-psi[i-1] )
}


plot(psi1,type = 'l')

psi3 = matrix(0,length(psi),2)
psi3[,2]=psi1
psi3[,1]=psi
k = rep(0,2)
for (i in 1: length(psi))
{
    for (j in (i):length(psi))
    {
        if (psi3[i,1]>psi3[j,1])
        {
            k = psi3[i,]
            psi3[i,]=psi3[j,]
            psi3[j,]=k
        }
    }
}
plot(psi3,type='l')


#plot(ts,psi,type = 'h')
#plot.spec.phase(atan2(Im(c),Re(c)))

#Фазовый спектр
#psi = atan2(Im(c),Re(c))
psi = matrix(0,length(c),2)
psi[,1] = 1:length(c)
psi[,2] = abs((1 - 1:length(c) *1i)/(1+(1:length(c))^2))
plot(psi,type = 'h')

#производная
psi1 = rep(0,length(psi))
for (i in 2:length(psi))
{
    psi1[i] = psi[i]-psi[i-1] 
}
plot(psi1,type = 'h')



psi2 = fft(c)
psi3 = atan2(Im(psi2),Re(psi2))
#производная
psi1 = rep(0,length(psi))
for (i in 2:length(psi))
{
    psi1[i] = psi3[i]-psi3[i-1] 
}
plot(psi1,type = 'h')






f = abs(c)

acq.freq <- length(c)                    # data acquisition (sample) frequency (Hz)
time     <- 4                      # measuring time interval (seconds)
ts       <- seq(0,time-1/acq.freq,1/acq.freq) # vector of sampling time-points (s) 
f.0 <- 1/time

dc.component <- 1
component.freqs <- c(3,7,10)        # frequency of signal components (Hz)
component.delay <- c(0,0,0)         # delay of signal components (radians)
component.strength <- c(1.5,.5,.75) # strength of signal components

f   <- function(t,w) { 
    dc.component + 
        sum( component.strength * sin(component.freqs*w*t + component.delay)) 
}

plot.fourier(f,f.0,ts=ts)



c1 = fft(revE,inverse = T)
c2 = fft(revE)
c3 = (c1[2:length(c1)]-c2[2:length(c2)])/2
#c3 = Im(c3)
psi = c3

psi1 = rep(0,length(psi))
for (i in 2:length(psi))
{
    psi1[i] = psi[i]-psi[i-1] 
}
plot(psi1[1:(length(psi1)/2)],type = 'l')


#-------------------
#Шарий1

#Предполагается, что сигнал на подасчу - одно слово
zcrm = function(signal,wlen)
{
    NumFrames = floor(length(signal)/wlen)
    zcrf = rep(0,NumFrames)
    for (j in 1:NumFrames)
    {
        for (i in ((j-1)*wlen+1):(j*wlen-1))
        {
            zcrf[j] = zcrf[j] + (sign(signal[i+1]) - sign(signal[i]))/2
        }
    }
    zcrf = zcrf/wlen
    return (zcrf)
}

for (j in 1:(length(words)/2))
{
signal = x1@left#[(words[2*j-1]*0.05*x1@samp.rate):(words[2*j]*0.05*x1@samp.rate)]
signalm = matrix(0,length(signal),2)
signalm[,1]= (1:length(signal))/x1@samp.rate
signalm[,2] = signal


E = Energy(signal,0.02*x1@samp.rate,0.02*x1@samp.rate)
a = zcr(signal,x1@samp.rate,wl=0.02*x1@samp.rate,plot = F)

M.ZCR = mean(a[,2])
SD.ZCR = sd(a[,2])*sd(a[,2])

M.E = mean(E)
SD.E = sd(E)*sd(E)
for (i in 1:nrow(a))
{
    if (a[i,2]>(M.ZCR+0.5*SD.ZCR))
    {
        a[i,2]=5000
    } else
    {
        if (E[i] < M.E - 0.6*SD.E)
        {
            a[i,2] = -2500
        } else
        {
            a[i,2] = 5000
        }
    }
    
}

plot(signalm,type = 'l')
#lines(words,rep(100000,length(words)),col = 'green',type = 'h')
#lines(a,type = 'h',col = 'red')
b = matrix(0,nrow(a),2)
for (i in 4:(nrow(a)-4))
{
    if (a[i,2]<0 && a[i-1,2]<0 && a[i-2,2]<0 && a[i-3,2]<0 && a[i+1,2]>0
        && (a[i+2]>0 | a[i+3]>0 | a[i+4]>0))
    {
        b[i,1] = a[i,1]
        b[i,2] = 30000 
    }
}
lines(b,type = 'h',col = 'green')
}
#Видимо условия 1) наличие небольшого количества гласных 2) остутсвтие согласных







#--------------------- еще один способ фрагментации
for (j in 1:(length(words)/2))
{
signal = x1@left[(words[j*2-1]*0.05*samprate):(words[j*2]*0.05*samprate)]
signal = signal/max(signal)

windowlen = 0.01*samprate
steplen = windowlen/2
NumFrames = floor((length(signal)-windowlen)/steplen) + 1
pos = 1
E = matrix(0,NumFrames,6)

#Энергия по каждому вейвлету
for (i in 1:(NumFrames-1))
{
    window = signal[pos:(pos+windowlen-1)]
    pos = pos + steplen
    for (k in 1:length(window))
    {
        window[k] = window[k] * (0.53836 - 0.46164 * cos (2 * pi * k / (length(window)-1) ))
    }
    dwtsignal = dwt (window, n.levels = 6)
    E[i,1] = sum((dwtsignal@W$W1)^2)
    E[i,2] = sum((dwtsignal@W$W2)^2)
    E[i,3] = sum((dwtsignal@W$W3)^2)
    E[i,4] = sum((dwtsignal@W$W4)^2)
    E[i,5] = sum((dwtsignal@W$W5)^2)
    E[i,6] = sum((dwtsignal@W$W6)^2)
}

#Сглаживаем
for (i in 1:NumFrames)
{
    for (l in 1:3) 
    {
        for (k in 1:floor(nrow(E)/3))
        {
            E[(3*k-2):(3*k),l] = max(E[(3*k-2):(3*k),l])
        }
    }
    for (l in 4:6) 
    {
        for (k in 1:floor(nrow(E)/5))
        {
            E[(5*k-4):(5*k),l] = max(E[(5*k-4):(5*k)],l)
        }
    }
}

#Производные
R = matrix(0,(NumFrames-1),6)
for (l in 1:6)
{
    R[,l] = diff(E[,l])/2
}

DIVopt = 0.05
EMin = 0.1/500

pause = matrix(1,6*NumFrames,2)
pos = 1
for (k in 1:1)
{
for (i in 3:(nrow(R)-2))
{
    if (#(abs(abs(R[i,k]) - E[i,k])<DIVopt)
        #&& (
        #    (E[i-1,k] - abs(R[i-1,k]) > 0.5 *  DIVopt) | (E[i+1,k] - abs(R[i+1,k]) >  0.5 * DIVopt)
        #    )
         #&& 
         (E[i,k] > EMin)#(mean(E[,k])-0.75*sd(E[,k])))
    )
    {
        pause[pos,1] = i
        pos = pos + 1
    }
}
}
pause[,1] = pause[,1]*(steplen+1)

plot(signal,type = 'l')
lines(pause,type = 'h',col = 'red')
lines(((1:length(E[,1])*(steplen+1))),E[,1]*500,type = 'l',col = 'green')

}



#r = kmeans(p,4)
#r$centers*(steplen)
