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
library(TTR)
setwd ("/Users/Kondor/Desktop/Music/FromC#")
x1 = readWave("9.wav")
x1 = updateWave(x1)

maxfreq = 512
sample_size = floor(maxfreq/2) #  

#floor((length(x1@left) - sample_size)/44100-3)))
print (Sys.time())
nframe = 17000 # до 13 минут
#magnitude = matrix(0,1,nframe*maxfreq)
magnitude = matrix(0,nframe,1)
for (j in (0):(nframe-1)) 
{
    sample1 = (j * sample_size)
    sample2 = (j + 1) * sample_size
    if (sample1 > 0 & sample2 < length(x1@left))
    {
        frame1 = x1@left[sample1:sample2]
        q = fft(frame1, 1)
        q = q/length(q)
        #q = q[300:length(q)]
        magnitude1 = (Re(q)*Re(q) + Im(q)*Im(q))^0.5 #амплитуда
        magnitude[j] = max(magnitude1)
    }
}
#-----------------------конец цикла
print (Sys.time())
magnitude = magnitude[magnitude != 0]
magnitud = magnitude[magnitude != 0]
plot(magnitude,type = 'h')

lines (EMA(magnitude,70), col = 'red')
lines (SMA(magnitude,70), col = 'green')


x1@samp.rate
 

#Нужно избавиться от тонких пиков


#выделяем время слога  !первый сигнал нужно будет отбросить, он просто регистрирует начало записи
locmax = 0 # переменная, отвечающая за локальный максимум
locmaxi = 0 # координата текущего locmax
for ( i in 21:(length(magnitude)-5))
    if ((3*magnitude[i-5]<magnitude[i] | 3*magnitude[i-6]<magnitude[i] | 
        3*magnitude[i-7]<magnitude[i] | 3*magnitude[i-8]<magnitude[i] |
        3*magnitude[i-9]<magnitude[i] | 3*magnitude[i-10]<magnitude[i]) &&  #Условие на максимум
        (3*magnitude[i+3] > magnitude[i])) #защита от одиночных случайных пиков
    {
        if (magnitude[i] > locmax)
        {
            magnitude1[i] = magnitude[i]
            locmax = magnitude[i]
            magnitude1[locmaxi] = 0
            locmaxi = i
        } else
        {
            magnitude1[i] = 0
        }
    }else #значит, пошел следующий слог
    {
        magnitude1[i] = 0  
        locmax = 0
        locmaxi = 0
    }

for (i in 11:(length(magnitude1)))
{
    if (magnitude1[i] != 0)
    {
        magnitude1[i] = max(magnitude1[i:(i-10)])
        magnitude1[(i-1):(i-10)] = 0
        if (magnitude1[i] < 100) #убираем шум
            magnitude1[i] = 0   
    }
}
    

plot(magnitude1,type = 'h')
steptime = matrix (0,length(magnitude1[magnitude1 != 0]),2)
a = 0
for (i in 10:length(magnitude1))
    if (magnitude1[i] != 0)
    {
        a = a + 1
        steptime[a,1] = i*sample_size /x1@samp.rate
        steptime[a,2] = magnitude[i]
        
    }

#Список время/амплитуда
steptime
