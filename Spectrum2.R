#-----------------------------------------------------------
#fft

SinglePi = pi
DoublePi = 2*pi

FFT_time = function (frame, direct)
{
    if (length(frame) == 2)
    { 
        spectrum = complex (0i,2) 
        spectrum[1] = frame[1] + frame[2] 
        spectrum[2] = frame[1] - frame[2]
        FFT_time = spectrum
    } 
    else{
        frameHalfSize = floor(length(frame)/2) # frame.Length/2
        frameFullSize = length(frame)
        
        frameOdd = complex (0i,frameHalfSize) 
        frameEven = complex (0i,frameHalfSize) 
        
        for ( i in  1:frameHalfSize)
        {
            frameOdd[i] = frame[2*i] 
            frameEven[i] = frame[2*i-1] 
        }
        
        spectrumOdd = FFT_time(frameOdd, direct) 
        spectrumEven = FFT_time(frameEven, direct) 
        
        if ( direct ==1){
            arg = -DoublePi/frameFullSize
        }else{
            arg = DoublePi/frameFullSize
        }
        
        omegaPowBase = cos(arg) + 1i * sin(arg) 
        omega = 1+0i 
        
        spectrum = complex (0i,frameFullSize) 
        
        for (j in 1:frameHalfSize)
        {
            spectrum[j] = spectrumEven[j] + omega*spectrumOdd[j]   #Внимание! комплексное произведение!
            spectrum[j + frameHalfSize] = spectrumEven[j] - omega*spectrumOdd[j] 
            omega = omega*omegaPowBase 
        }
        
        FFT_time = spectrum 
    }}


#-----------------------------------------------------------

#Вспомогательная функция
Align = function (angle, period)
{
    angle = (angle +pi) %% period - pi
    Align = angle
    #qpd = as.integer (angle/period);
    #if (qpd >= 0) {
    #    qpd =qpd + qpd %% 2
    #} else{
    #    qpd = qpd - qpd %% 2
    #}
    #Align = angle - period*qpd
}


#Функция уточнения сигнала на основе сдвига
SpectrumJoin = function (spectrum1, spectrum2, shiftsPerFrame, sampleRate)
{
    frameSize = length(spectrum1) #длина фрейма
    frameTime = frameSize/sampleRate #время в секундах
    shiftTime = frameTime/shiftsPerFrame #время сдвига в секундах
    binToFrequancy = sampleRate/frameSize 
    dictionary =  matrix(0,frameSize,2)
    
    magnitude1 = (Re(spectrum1)*Re(spectrum1) + Im(spectrum1)*Im(spectrum1))^0.5 #амплитуда
    phase1 = atan2(Im(spectrum1),Re(spectrum1)) #Арктангенс (y,x) #фаза
    magnitude2 = (Re(spectrum2)*Re(spectrum2) + Im(spectrum2)*Im(spectrum2))^0.5 #амплитуда
    phase2 = atan2(Im(spectrum2),Re(spectrum2)) #Арктангенс (y,x) #фаза
    magnitude0 = magnitude1 + magnitude2 
    
    for ( bin in  1:frameSize)
    {
        omegaExpected = DoublePi*(bin*binToFrequancy) # ω=2πf
        omegaActual = (phase2[bin] - phase1[bin])/shiftTime # ω=∂φ/∂t  
        #print(omegaActual)
        #print(omegaExpected)
        omegaDelta = Align(omegaActual - omegaExpected, DoublePi) # Δω=(∂ω + π)%2π - π  
        binDelta = omegaDelta/(DoublePi*binToFrequancy) 
        frequancyActual = (bin + binDelta)*binToFrequancy 
        #print(frequancyActual)
        dictionary[bin,1] = frequancyActual
        dictionary[bin,2] = magnitude0[bin]*(0.5 + abs(binDelta))
    }
    SpectrumJoin = dictionary
}

#----------------------------------------------



library(tuneR)
setwd ("/Users/Kondor/Desktop/Music/FromC#")
x1 = readWave("mysong4.wav")
x1 = updateWave(x1)

#если частота сигнала лежит где-то между границами шага ближе к середине 
#выйходит пик со «срезанной» вершиной и нам затруднительно сказать, что же там за частота.
# Соседние могут склеиться
#Решается с помощью аппроксимаций (как в NaiveBeat), но это малоточно
#Данный метод уточнения частоты сигнала, основан на вычислении задержки фаз у спектров двух кадров,
#наложенных друг на друга, но немного сдвинутых во времени.

ShiftsPerFrame = 16 # мы сдвигаем frame на 1/16
# первый frame
maxfreq = 4096
sample_size = floor(maxfreq/2) #  

#floor((length(x1@left) - sample_size)/44100-3)))
qwerty <- vector()
for (j in (-250):(250))
{
sample1 = ((length(x1@left)/2) + j * sample_size)
sample2 = ((length(x1@left)/2) + (j + 1) * sample_size)
#sample1 = floor((length(x1@left)/2) - sample_size)
#sample2 = floor((length(x1@left)/2) + sample_size)
frame1 = x1@left[sample1:sample2]

#применяем окно Hann -дифракционные эффкты убираем (сейчас - Гаусс)

for (i in 1:length(frame1))
{
    a = (length(frame1) - 1)/2;
    t = (i - a)/(0.5*a);
    t = t*t
    exp(-t/2);
    frame1[i] = frame1[i]*exp(-t/2);
}
spectrum1 = FFT_time(frame1, 1);

#Второй фрейм
#sample1 = floor((length(x1@left)/2) - sample_size*(ShiftsPerFrame-1)/ShiftsPerFrame)
#sample2 = floor((length(x1@left)/2) + sample_size*(ShiftsPerFrame+1)/ShiftsPerFrame)

sample1 = ((length(x1@left)/2) + j * sample_size)
sample2 = ((length(x1@left)/2) + (j + 1) * sample_size)
frame2 = x1@left[sample1:sample2]
#применяем окно Hann -дифракционные эффкты убираем (сейчас - Гаусс)

for (i in 1:length(frame2))
{
    a = (length(frame1) - 1)/2;
    t = (i - a)/(0.5*a);
    t = t*t
    exp(-t/2);
    frame2[i] = frame2[i]*exp(-t/2);
}
spectrum2 = FFT_time(frame2, 1);

#Нормируем
spectrum1 = spectrum1/length(spectrum1)
spectrum2 = spectrum2/length(spectrum2)

#на вход - спекты двух фреймов, на сколько сдвиг, частота дискретизации (HZ)
#В результате мы получим словарь частота-амплитуда, 
#где значения частот будут довольно близки к реальным.
SampleRate = x1@samp.rate
FinalSpectrum = SpectrumJoin(spectrum1, spectrum2, ShiftsPerFrame, SampleRate)


# то, что дальше, не написано до конца
#--------------------------------------


NOTE<-vector()
NOTE[1] = 0
LEN = vector()
LEN[1] = 1
a=1 #номер ноты
b=1 #Текущая длина ноты
for ( i in 2:length(FinalSpectrum[,1]))
{
    if (FinalSpectrum[i,2] > 50000 )#Проверка на значимость
    {
        if ((floor((log2(FinalSpectrum[i,2]) - log2(27.5))*12)) != (NOTE[a])) # проверка на новизну
        {
            LEN[a] = b
            a=a+1
            NOTE[a]=floor((log2(FinalSpectrum[i,1]) - log2(27.5))*12)
            b = 1
        } else {
            b=b+1
        }
        
    }else{
        b=b+1
    }
}

#plot(FinalSpectrum[,2],type ='l')
#plot(FinalSpectrum[50:500,2],type ='l')


z = 0
for (i in 1:length(FinalSpectrum[,1]))
{
    if (FinalSpectrum[i,2] == max(FinalSpectrum[,2]))
        z = FinalSpectrum[i,1]
}
qwerty[j] = floor((log2(z) - log2(27.5))*12)
}
#-----------------------конец цикла

NOTE<-vector()
NOTE[1] = qwerty[1]
LEN = vector()
LEN[1] = 1
a=1 #номер ноты
b=1 #Текущая длина ноты
for ( i in 2:length(qwerty))
{

    if (qwerty[i] != (NOTE[a])) # проверка на новизну
        {
            LEN[a] = b
            a=a+1
            NOTE[a]=qwerty[i]
            b = 1
        } else {
            b=b+1
        }
}
LEN[a] = b

b = data.frame(NOTE, LEN, T, T)
names(b) = c("note", "lenght","punctate","slur")

#qw1 = vector()
#for ( i in 0: (length(qwerty)-1))
#{
#    for (j in 1:(2*4096))
#    {
#       qw1[i*2*4096 + j] = qwerty[i+1]
#    }
#}
#qw2 = 27.5 * 2 ^ (qw1/12)
#writeWave(Wave(as.integer(qw2),bit = 24,samp.rate =44100), filename = "MyWave.wav")
#lilyinput(b, file = "Rsong.ly", midi =T)

