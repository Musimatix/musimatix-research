library(tuneR)
library(seewave)
setwd ("/Users/Kondor/Git/musimatix-research")
epsilon = 0.0000001
setwd ("/Users/Kondor/Desktop/Music/FromC#/korpus")
x1 = readWave("Tolya-Rodina-kuplet.wav")
x1 = updateWave(x1)
signal = x1@left


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
    Energyf = Energyf/sum(Energyf)
    return (Energyf)
}

windowlen = 0.05 #(mc)
steplen = windowlen # lenght of step

E = Energy(signal,windowlen*x1@samp.rate,steplen*x1@samp.rate)*100 #*100
Ewindowlen = windowlen*x1@samp.rate #lenght of window, which we use in Energy detection

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
    #if ((begin != 0) && (i-begin < 10) && (E[i]<0.1)) #too short word
    #    begin = 0
}
words = words[words != 0]
wordsh = matrix(-20000,length(words),2)
wordsh[,1] = words*windowlen*x1@samp.rate
words

#---------------------------------------------------

#generate points of local min and max of signal (fhi)
minmax = function (fhi)
{
    minimax = matrix(1,(length(fh)/5))
    steps = 1
    #распознавание последнего пика так как после может не быть спада спада
    fh = rep (0,(length(fhi)+15))
    for (i in 1:(length(fhi)))
        fh[i+5] = fhi[i]    #                                whats up??????????
    # Условие на макс
    for (i in 6:(length(fh)-5))
    {
        if (fh[i]>fh[i-1] &&fh[i]>fh[i-2] &&fh[i]>fh[i-3] &&fh[i]>fh[i-4] &&fh[i]>fh[i-5] &&
            fh[i]>fh[i+1] &&fh[i]>fh[i+2] &&fh[i]>fh[i+3] &&fh[i]>fh[i+4] &&fh[i]>fh[i+5] )
        {
            minimax[steps] = i-5
            steps = steps + 1
        }
    }
    # Условие на мин
    for (i in 6:(length(fh)-5))
    {
        if (fh[i]<= fh[i-1] &&fh[i] <= fh[i-2] &&fh[i] <= fh[i-3] &&fh[i] <= fh[i-4] &&fh[i] <= fh[i-5] &&
            fh[i] <= fh[i+1] &&fh[i] <= fh[i+2] &&fh[i] <= fh[i+3] &&fh[i] <= fh[i+4] &&fh[i] <= fh[i+5] 
            && (fh[i+1]>epsilon) #последнее условие из-за большого количества нулей
        )
        {
            minimax[steps] = i-5
            steps = steps + 1
        }
    }
    minimax = minimax[minimax != 1]
    return (minimax)
}

#функция-условие выбора слогов
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
    
    minmaxi = rep (1,(length(minmaxii)+(2*length(minmaxii))))
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
    #бегаем по значимым точкам массива с краями 
    
    
    #for (i in ((floor(length(minmaxii)))+1):((floor(length(minmaxi)-length(minmaxii)))))
    #{
    #    par = 0
    #    
    #    for (l in (-floor(length(minmaxii))):(-1))
    #    {
    #        for (ll in 1:(floor(length(minmaxii)))){
    #            
    #            if ((fh[minmaxi[i]] > locmax) &&
    #                (max(fh[minmaxi[(i+l):(i-1)]]) < fh[minmaxi[i]]) && (max(fh[minmaxi[(i+1):(i+ll)]]) < fh[minmaxi[i]] ) 
    #                && (min(fh[minmaxi[(i+l):(i-1)]]) < locmax) && ( min(fh[minmaxi[(i+1):(i+ll)]])< locmax) 
    #                && (abs(min(fh[minmaxi[(i+l):(i-1)]]) - fh[minmaxi[i]]) > locmax) && ( abs(min(fh[minmaxi[(i+1):(i+ll)]]) - fh[minmaxi[i]]) > locmax) 
    #                (fh[minmaxi[(i+l)]] < locmax/2) && (fh[minmaxi[(i+ll)]] < locmax/2)
    #                 && ((fh[minmaxi[i+l]] - fh[minmaxi[i]])/(2*(minmaxi[i]-minmaxi[i+l])) < -0.5)
    #            )
    #            {
    #                par = 1
    #            }
    #        }}
    #    if (par == 1)
    #    {
    #        convexi[steps] = minmaxi[i]
    #        steps = steps +1
    #    } 
    #}
    locmax = fh[1]
    locmaxi = 1
    locmin = max(fh)
    locmini = 1
    flag = 1 # 1 - последняя точка - макс
    for (i in ((floor(length(minmaxii)))+1):((floor(length(minmaxi)-length(minmaxii)))))
    {
        #если предыдущая точка была макс, и у нас макс(лучше)
        if (fh[minmaxi[i]]>locmax && flag == 1)
        {
            #new max & delete old max
            convexi[steps] = minmaxi[i]
            locmax = minmaxi[i]
            locmaxi = i
        }
        #если предыдущая точка была мин и у нас макс (плюс проверка на два очень близких пика)
        if ((fh[minmaxi[i]]-fh[minmaxi[i-1]])>0.08 && flag == 0 )
        {
            steps = steps + 1
            convexi[steps] = minmaxi[i]
            locmax = fh[minmaxi[i]]
            flag = 1
        }
        #если предыдущая точка была макс и у нас мин
        if (flag == 1 && (fh[minmaxi[i-1]] - fh[minmaxi[i]] > 0.1))
        {
            flag = 0
            locmax = 0
        }
        #если предыдущая точка была мин и у нас мин
        if (flag == 0 && (minmaxi[i-1] > minmaxi[i]))
        {
            convexi[steps] = minmaxi[i]
        }
    }
    return(convexi[convexi !=0 ])
}

#-------------------------
# Другие окна!
wlen = 0.01*x1@samp.rate
step = wlen/2
pause = matrix(30000,1000,2) #подумать!!!!!!!!!!!!!!!!!!!!!!!!
t = 1 #первый пропускаем

y1 = Sys.time()
#Рассматриваем только участки с большой энергией, что бы не было шумов
ME = mean(E)/1.2
for ( i in 1:(length(signal)-Ewindowlen))
{
    #print (E[floor(i/Ewindowlen)+1])
    if (E[floor(i/Ewindowlen)+1]< (ME)) #!!!!!!!!!!!!!!!!!!!  подобрать коэфф!!
    {
        signal[i] = 0
    }
}

y2 = Sys.time()
#Бегаем по каждому слову
#for (j in 1:(length(words)/2))
for (j in 4:4)
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
            if (st > 300 && st < 700 )
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
    fh = log(fh[20:(length(fh)-20)]+1)
    fh = fh/max(fh)
    a = minmax(fh)
    a = sort(a)
    b = convex(a,fh,0.2)
    if (length( b) != 0){
        b = b+20
        for (i in 1:length(b))
        {
            if (((words[2*j-1]*windowlen*x1@samp.rate+(b[i]*step + step) - pause[t-1])>0.15*x1@samp.rate))
            {
                pause[t,1] = words[2*j-1]*windowlen*x1@samp.rate+(b[i]*step + step)
                t= t+1
            }
        }}
    
}

y3 = Sys.time()
#---------------------------------
#Проверка
fhii = log(fh[20:(length(fh)-20)]+1)
fhii = fhii/max(fhii)
fhii = fh
aa = minmax(fhii)
aa = sort(aa)
plot(fhii,type = 'h')

bb = convex(aa,fhii,0.2)
cc = matrix(2, length(bb),2)
cc[,1]=bb
cc[,2] = 8
dd = matrix(2, length(aa),2)
dd[,1]=aa
dd[,2] = 0.5


points (dd,col = 'red')
points (cc, type = 'h',col = 'red')

plot(signal,type = 'l')
points(wordsh,col = 'green',type = 'h')
points(pause,type = 'h',col = 'red')

y1
y2
y3
