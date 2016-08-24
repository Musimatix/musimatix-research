library(tuneR)
library(seewave)
epsilon = 0.0000001
setwd ("/Users/Kondor/Desktop/Music/FromC#/korpus")
x1 = readWave("Tolya-Rodina-pripev.wav")
x1 = updateWave(x1)

q0 = Sys.time()



#производная
derivety = function (fun)
{
    derivetyf = rep (0,length(fun))
    derivetyf[1] = 0
    for (i in 2:length(fun))
        derivetyf[i] = (fun[i]-fun[i-1])/2
    return (derivetyf)
}


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



windowlen = 0.05
steplen = windowlen

E = Energy(x1@left,windowlen*x1@samp.rate,steplen*x1@samp.rate)*100
Ewindowlen = windowlen*x1@samp.rate
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

plot(E,type = 'l')
points(beatmin,col = 'green',type = 'h')
q1 = Sys.time()

#--------------------------
words = matrix(0,countlocmax)
begin = 0
step = 1
#Выделение слов
for ( i in 1:length(E))
{
    if ((begin == 0) && (E[i]>0.1)) #начало возможного слова
    {
        if (step != 1){
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
wordf = matrix(-20000,length(words),2)
wordf[,1] = words*windowlen*x1@samp.rate
points(wordf,type = 'h')
words
#Анализ слов по отдельности
q2 = Sys.time()
#оставим только сигнал слова
signal = x1@left


# number of frame, fh - magnitude ~500hz
fhbool1 = function(i,fh)
{
    poss = 0
    
        if (fh[i] < (1*mean(fh)) && fh[i-1] < (mean(fh)) && fh[i-2] < (mean(fh)) 
            && fh[i-3] < (mean(fh)) && fh[i-4] < (mean(fh)) && fh[i-5] < (mean(fh))
            && fh[i-6] < (mean(fh)) && fh[i-7] < (mean(fh)) && fh[i-8] < (mean(fh))
            && fh[i-9] < (mean(fh)) && fh[i-10] < (mean(fh)) && fh[i+1] > mean(fh))
            poss = poss + 1
        
        if (poss > 0)
        {
            #print (i)
            return (TRUE)
        } else
        {
            return (FALSE)
        }
}

minmax = function (fhi)
{
    minimax = matrix(1,(length(fh)/10))
    steps = 1
    #распознавание последнего пика так как после может не быть спада спада
    fh = rep (0,(length(fhi)+15))
    for (i in 1:(length(fhi)))
         fh[i] = fhi[i]
    # Условие на макс
    print (fh)
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
            && (fh[i+1]>epsilon)
        )
        {
            minimax[steps] = i
            steps = steps + 1
        }
    }
    #minimax = minimax + 6
    print (minimax)
    return (minimax)
}

convex = function (minmaxii,fh,d)
{
    convexi = matrix (0,(length(minmaxii)/2))
    
    #min of minmax
    locmin = max(fh[minmaxii],na.rm = T)
    locmini = 1

    for (i in 1:(length(minmaxii)))
    {
       
        if (fh[minmaxii[i]]< locmin)
        {
            locmin = fh[minmaxii[i]]
            locmini = i
        }
    }
    
    minmaxi = rep (minmaxii[locmini],(length(minmaxii)+(2*length(minmaxii))))
    for (i in 1:length(minmaxii))
    {
        minmaxi[i+(length(minmaxii))] = minmaxii[i]
    }
    
    #max of minmax что бы избавиться от самого большого пика в оценке
    locmax = min(minmaxii)
    locmaxi = 1
    for (i in 1:length(minmaxii))
    {
        if (fh[minmaxii[i]]> locmax)
        {
            locmax = fh[minmaxii[i]]
            locmaxi = i
        }
    }
    steps = 1
    locmax = max(fh[minmaxii[-locmaxi]])/3
    #бегаем по значимым точкам массива с краями 
    for (i in ((floor(length(minmaxii)))+1):((floor(length(minmaxi)-length(minmaxii)))))
    {
        #if (
            #(fh[minmaxi[i]] - fh[minmaxi[i+1]] > d) && ((fh[minmaxi[i]] - fh[minmaxi[i-1]] > d)) |
            #(fh[minmaxi[i]] - fh[minmaxi[i+2]] > d) && ((fh[minmaxi[i]] - fh[minmaxi[i-1]] > d)) |
            #(fh[minmaxi[i]] - fh[minmaxi[i+3]] > d) && ((fh[minmaxi[i]] - fh[minmaxi[i-1]] > d)) |
            #(fh[minmaxi[i]] - fh[minmaxi[i+1]] > d) && ((fh[minmaxi[i]] - fh[minmaxi[i-2]] > d)) |
            #(fh[minmaxi[i]] - fh[minmaxi[i+2]] > d) && ((fh[minmaxi[i]] - fh[minmaxi[i-2]] > d)) |
            #(fh[minmaxi[i]] - fh[minmaxi[i+3]] > d) && ((fh[minmaxi[i]] - fh[minmaxi[i-2]] > d)) |
            #(fh[minmaxi[i]] - fh[minmaxi[i+1]] > d) && ((fh[minmaxi[i]] - fh[minmaxi[i-3]] > d)) |
            #(fh[minmaxi[i]] - fh[minmaxi[i+2]] > d) && ((fh[minmaxi[i]] - fh[minmaxi[i-3]] > d)) |
            #(fh[minmaxi[i]] - fh[minmaxi[i+3]] > d) && ((fh[minmaxi[i]] - fh[minmaxi[i-3]] > d)) |
            #(fh[minmaxi[i]] - fh[minmaxi[i+1]] > 0.5*d) && ((fh[minmaxi[i]] - fh[minmaxi[i-1]] > 2*d)) |
            #(fh[minmaxi[i]] - fh[minmaxi[i+1]] > 0.5*d) && ((fh[minmaxi[i]] - fh[minmaxi[i-2]] > 2*d)) |
            #(fh[minmaxi[i]] - fh[minmaxi[i+1]] > 0.5*d) && ((fh[minmaxi[i]] - fh[minmaxi[i-3]] > 2*d)) 
            
            #((fh[minmaxi[i]] > locmax) && (fh[minmaxi[i-1]] < locmax) && (fh[minmaxi[i+1]] < locmax)) |
            #((fh[minmaxi[i]] > locmax) && (fh[minmaxi[i-1]] < locmax) && (fh[minmaxi[i+2]] < locmax) && (fh[minmaxi[i+1]] < fh[minmaxi[i]] )) |
            #((fh[minmaxi[i]] > locmax) && (fh[minmaxi[i-1]] < locmax) && (fh[minmaxi[i+3]] < locmax) && (fh[minmaxi[i+1]] < fh[minmaxi[i]] ) && (fh[minmaxi[i+2]] < fh[minmaxi[i]] )) |
            #((fh[minmaxi[i]] > locmax) && (fh[minmaxi[i-2]] < locmax) && (fh[minmaxi[i+1]] < locmax) && (fh[minmaxi[i-1]] < fh[minmaxi[i]] )) |
            #((fh[minmaxi[i]] > locmax) && (fh[minmaxi[i-3]] < locmax) && (fh[minmaxi[i+1]] < locmax) && (fh[minmaxi[i-1]] < fh[minmaxi[i]] ) && (fh[minmaxi[i-2]] < fh[minmaxi[i]] )) |
            #((fh[minmaxi[i]] > locmax) && (fh[minmaxi[i-2]] < locmax) && (fh[minmaxi[i+2]] < locmax) && (fh[minmaxi[i+1]] < fh[minmaxi[i]] ) && (fh[minmaxi[i-1]] < fh[minmaxi[i]] )) |
            #((fh[minmaxi[i]] > locmax) && (fh[minmaxi[i-2]] < locmax) && (fh[minmaxi[i+3]] < locmax) && (fh[minmaxi[i+1]] < fh[minmaxi[i]] ) && (fh[minmaxi[i+2]] < fh[minmaxi[i]] ) &&(fh[minmaxi[i-1]] < fh[minmaxi[i]] )) |
            #((fh[minmaxi[i]] > locmax) && (fh[minmaxi[i-3]] < locmax) && (fh[minmaxi[i+2]] < locmax) && (fh[minmaxi[i+1]] < fh[minmaxi[i]] ) && (fh[minmaxi[i-1]] < fh[minmaxi[i]] ) &&(fh[minmaxi[i-2]] < fh[minmaxi[i]] )) |
            #((fh[minmaxi[i]] > locmax) && (fh[minmaxi[i-3]] < locmax) && (fh[minmaxi[i+3]] < locmax) && (fh[minmaxi[i+1]] < fh[minmaxi[i]] ) && (fh[minmaxi[i+2]] < fh[minmaxi[i]] ) &&(fh[minmaxi[i-1]] < fh[minmaxi[i]] ) && (fh[minmaxi[i-2]] < fh[minmaxi[i]] )) |
            
            #((fh[minmaxi[i]] > locmax) && (fh[minmaxi[i-1]] < locmax) && (fh[minmaxi[i+3]] < locmax) && (fh[minmaxi[i+1]] < fh[minmaxi[i]] ) && (fh[minmaxi[i+2]] < fh[minmaxi[i]] )) |
            #((fh[minmaxi[i]] > locmax) && (fh[minmaxi[i-1]] < locmax) && (fh[minmaxi[i+4]] < locmax) && (fh[minmaxi[i+1]] < fh[minmaxi[i]] ) && (fh[minmaxi[i+2]] < fh[minmaxi[i]] )) |
            #((fh[minmaxi[i]] > locmax) && (fh[minmaxi[i-1]] < locmax) && (fh[minmaxi[i+5]] < locmax) && (fh[minmaxi[i+1]] < fh[minmaxi[i]] ) && (fh[minmaxi[i+2]] < fh[minmaxi[i]] )) |
            #((fh[minmaxi[i]] > locmax) && (fh[minmaxi[i-1]] < locmax) && (fh[minmaxi[i+6]] < locmax) && (fh[minmaxi[i+1]] < fh[minmaxi[i]] ) && (fh[minmaxi[i+2]] < fh[minmaxi[i]] )) 
            
        #)
        par = 0
        
        for (l in (-floor(length(minmaxii))):(-1))
        {
            for (ll in 1:(floor(length(minmaxii)))){
                
                if ((fh[minmaxi[i]] > locmax) &&
                    (max(fh[minmaxi[(i+l):(i-1)]]) < fh[minmaxi[i]]) && (max(fh[minmaxi[(i+1):(i+ll)]]) < fh[minmaxi[i]] ) &&
                    (min(fh[minmaxi[(i+l):(i-1)]]) < locmax) && ( min(fh[minmaxi[(i+1):(i+ll)]])< locmax) &&
                    (abs(min(fh[minmaxi[(i+l):(i-1)]]) - fh[minmaxi[i]]) > locmax) && ( abs(min(fh[minmaxi[(i+1):(i+ll)]]) - fh[minmaxi[i]]) > locmax) &&
                    (fh[minmaxi[(i+l)]] < locmax/2) && (fh[minmaxi[(i+ll)]] < locmax/2)
                )
                {
                    par = 1
                }
        }}
        if (par == 1)
        {
            convexi[steps] = minmaxi[i]
            steps = steps +1
        } 
    }
    return(convexi[convexi !=0 ])
}

q3 = Sys.time()
wlen = 0.01*x1@samp.rate
step = wlen/2
pause = matrix(30000,1000,2) #подумать!!!!!!!!!!!!!!!!!!!!!!!!
t = 1

#Рассматриваем только участки с большой энергией, что бы не было шумов
for ( i in 1:(length(signal)-Ewindowlen))
{
    #print (E[floor(i/Ewindowlen)+1])
    if (E[floor(i/Ewindowlen)+1]< (max(E)/10)) #!!!!!!!!!!!!!!!!!!!  подобрать коэфф!!
    {
        signal[i] = 0
    }
}

#for (j in 1:(length(words)/2))
for (j in 3:3)
{
    #начала слов
    pause[t,1] = (words[2*j-1]*windowlen*x1@samp.rate)
    t = t+1
    #нормируем
    wsignal = signal[(words[2*j-1]*windowlen*x1@samp.rate):(words[2*j]*windowlen*x1@samp.rate)]
    wsignal = wsignal/max(wsignal)
    len = length(wsignal)
    pos = 1
    #Число фреймов
    NumFrames = floor((len-wlen)/step) + 1
    fh = matrix(0,NumFrames) #пятисотая
    for (i in 1:NumFrames)
    {
        window = wsignal[pos:(pos+wlen-1)]
        pos = pos + step
        x = fft(window)
        bin = x1@samp.rate/length(x)
        st = 0
        for (k in 1:(length(x)/2))
        {
            st = st + bin
            if (st > 500 && st < 2000 )
            {
                fh[i] = abs(x[k])
            }
        }
    }
    #просто удобная нормировка
    #fh = log10 (fh+1)
    fh1 = fh
    for (i in 8:(length(fh)-13))
        fh1[i]=1/21 * sum(fh[(i-7):(i+13)])
    fh=fh1
    a = minmax(fh[20:(length(fh)-20)])
    a = sort(a)
    b = convex(a,fh[20:(length(fh)-20)],0.2)
    if (length( b) != 0){
    print (b+20)
    b = b+20
    #plot(fh,type = 'h')
    for (i in 1:length(b))
    {
        
        if (((words[2*j-1]*windowlen*x1@samp.rate+(b[i]*step + step) - pause[t-1])>0.15*x1@samp.rate))
        {
            pause[t,1] = words[2*j-1]*windowlen*x1@samp.rate+(b[i]*step + step)
            t= t+1
        }
    }}
    
}

q4 = Sys.time()
q0
q1
q2
q3
q4

head(pause,10)
#plot(signal,type = 'l')
#points(wordf,col = 'green',type = 'h')
#points(pause,type = 'h',col = 'red')


aa = minmax(fh[20:(length(fh)-20)])
aa = sort(aa)


plot(fh[20:(length(fh)-20)],type = 'h')

bb = convex(aa,fh[20:(length(fh)-20)],0.2)
cc = matrix(2, length(bb),2)

cc[,1]=bb
cc[,2] = 8

dd = matrix(2, length(aa),2)
dd[,1]=aa
dd[,2] = 0.5


points (dd,col = 'red')
points (cc, type = 'h',col = 'red')



#for (j in 1:(length(words)/2))
# - это номер слова
for (j in 2:2)
{
    #нормируем
    print (j)
    wsignal = signal[(words[2*j-1]*windowlen*x1@samp.rate):(words[2*j]*windowlen*x1@samp.rate)]
    wsignal = wsignal/max(wsignal)
    len = length(wsignal)
    pos = 1
    #Число фреймов
    NumFrames = floor((len-wlen)/step) + 1
    fh = matrix(0,NumFrames) #низкние частоты
    for (i in 1:NumFrames)
    {
        window = wsignal[pos:(pos+wlen-1)]
        for (k in 1:wlen)
        {
            window[k] = window[k] * (0.53836 - 0.46164 * cos (2 * pi * k / (wlen - 1)))
        }
        
        pos = pos + step
        x = fft(window)
        bin = x1@samp.rate/length(x)
        st = 0
        for (k in 1:(length(x)/2))
        {
            st = st + bin
            if (st > 60 && st < 300 )
            {
                fh[i] = abs(x[k])
            }
        }
    }
    
    #просто удобная нормировка
    fh1 = fh
    for (i in 8:(length(fh)-13))
        fh1[i]=1/21 * sum(fh[(i-7):(i+13)])
    fh=fh1
    #plot (fh,type = 'h')
    paus = matrix (0,length(pause[(pause[,1] != 0),1]))
    paus = pause[(pause[,1] != 0),1]# zeros out
    pauswin = matrix (0,length(pause[(pause[,1] != 0),1]))
    for ( i in 1:(length(paus)))
    {
        if (paus[i] <= (words[2*j]*windowlen*x1@samp.rate)  && paus[i] >= (words[2*j-1])*windowlen*x1@samp.rate)
        {
            pauswin[i] = paus[i]
        }
    }
    pauswin = pauswin[pauswin !=0] # zeros out
    
    #fundamential freq
    #1 - position in whole signal, 2 - max f0,3 - mean f0, 4 - nuber of up frames, 5 - number of down frames
    #6 - 4/5 7 - 3/4
    f0 = matrix(0,length(pauswin),7)
    
    #for every syll (except the last)
    if (length(pauswin) > 2)
    {
    for ( i in 1:(length(pauswin)-1))
    {
        syll = signal[pauswin[i]:pauswin[i+1]]
        x = abs(fft((syll))[1:(length(syll)/2)])
        #plot(x, type = 'l')
        locmax = 0
        locmaxi = 1
        for (k in 1:length(x))
        {
            if (x[k] > locmax)
            {
                locmax = x[k]
                locmaxi = k
            }
        }
        #координата слога
        f0[i,1] = pauswin[i]
        # max значение f0 на всем слоге
        f0[i,2] = locmax
        NumFrames = floor(length(syll)/0.01/x1@samp.rate)
        #mean (f0) in syll
        mf0 = matrix(0,1:NumFrames)
        poss = 1
        for (k in 1:NumFrames)
        {
            window = syll[(poss):(poss+0.01*x1@samp.rate-1)]
            poss = poss+1
            x = abs(fft((window))[1:(length(window)/2)])
            #plot(x, type = 'l')
            locmax = 0
            for (l in 1:length(x))
            {
                if (x[l] > locmax)
                {
                    locmax = x[l]
                    locmaxi = l
                }
            }
            mf0[k] = locmax
        }
        f0[i,3] = log(mean(mf0)) #log of mean f0 frames in syll
        #number of ups
        ups = 0
        for ( k in 1:length(mf0))
        {
            if (derivety(mf0)[k] > 0)
                ups = ups + 1
        }
        f0[i,4] = ups
        
        #number of downs
        downs = 0
        for ( k in 1:length(mf0))
        {
            if (derivety(mf0)[k] < 0)
                downs = downs + 1
        }
        f0[i,5] = downs
        f0[i,6] = f0[i,4]/f0[i,5]
        f0[i,7] = f0[i,3]/f0[i,4]
        
        #/ (f0[i,1] - f0[i+1,1]) * 2 # psevdo slope
    }
        #the last syll
                i = length(pauswin)
                #тут слог от своего начала до конца слова
                syll = signal[pauswin[i]:(words[2*j]*windowlen*x1@samp.rate)] 
                x = abs(fft((syll))[1:(length(syll)/2)])
                #plot(x, type = 'l')
                locmax = 0
                locmaxi = 1
                for (k in 1:length(x))
                {
                    if (x[k] > locmax)
                    {
                        locmax = x[k]
                        locmaxi = k
                    }
                }
                #координата слога
                f0[i,1] = pauswin[i]
                # max значение f0 на всем слоге
                f0[i,2] = locmax
                NumFrames = floor(length(syll)/0.01/x1@samp.rate)
                #mean (f0) in syll
                mf0 = matrix(0,1:NumFrames)
                poss = 1
                for (k in 1:NumFrames)
                {
                    window = syll[(poss):(poss+0.01*x1@samp.rate-1)]
                    poss = poss+1
                    x = abs(fft((window))[1:(length(window)/2)])
                    #plot(x, type = 'l')
                    locmax = 0
                    for (l in 1:length(x))
                    {
                        if (x[l] > locmax)
                        {
                            locmax = x[l]
                            locmaxi = l
                        }
                    }
                    mf0[k] = locmax
                }
                f0[i,3] = log(mean(mf0)) #log of mean f0 frames in syll
                #number of ups
                ups = 0
                for ( k in 1:length(mf0))
                {
                    if (derivety(mf0)[k] > 0)
                        ups = ups + 1
                }
                f0[i,4] = ups
                
                #number of downs
                downs = 0
                for ( k in 1:length(mf0))
                {
                    if (derivety(mf0)[k] < 0)
                        downs = downs + 1
                }
                f0[i,5] = downs
                f0[i,6] = f0[i,4]/f0[i,5]
                f0[i,7] = f0[i,3]/f0[i,4]
                #/ (f0[i,1] - f0[i+1,1]) * 2 # psevdo slope
        }  
    
}
pauswin
pauswini = matrix(15000, length(pauswin),2)
pauswini[,1] = pauswin
f0
plot(signal,type = 'l')
points(wordf,col = 'green',type = 'h')
points(pause,type = 'h',col = 'red')
points(pauswini,col = 'green')


#Делаем предположение, что в слове не более 2ух ударных, иначе просто порежем длинные слова
