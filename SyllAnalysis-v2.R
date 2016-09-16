library(tuneR)
library(seewave)
library(TTR)
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


meanEnergy = function (signal, wlen, step)
{
    EN = Energy(signal, wlen, step)
    meanEnergyf = SMA(EN,50)
    meanEnergyf[1:50] = rep(0,50)
    return (meanEnergyf)
}

#meanEnergyq = meanEnergy(signal,windowlen*x1@samp.rate,steplen*x1@samp.rate)
#for (i in 1:length(E))
#{
#    if ((EN[i]==meanEnergyq[i]) && EN[i] != 0 )
#    {
#        print(i)
#    }
#}
    
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
    minimax = matrix(1,(length(fhi)/5))
    steps = 1
    #распознавание последнего пика так как после может не быть спада спада
    fh = rep (0,(length(fhi)+15))
    for (i in 1:(length(fhi)))
        fh[i+5] = fhi[i]    #                                whats up??????????
    # Условие на макс
    for (i in 6:(length(fh)-5))
    {
        if ((fh[i]>=fh[i-1]) && (fh[i]>=fh[i-2]) && (fh[i]>=fh[i-3]) && (fh[i]>=fh[i-4]) && (fh[i]>=fh[i-5]) &&
            (fh[i]>=fh[i+1]) && (fh[i]>=fh[i+2]) && (fh[i]>=fh[i+3]) && (fh[i]>=fh[i+4]) && (fh[i]>=fh[i+5]) 
            && (abs(fh[i] - fh[i-1]) > epsilon)) #последнее условие из-за большого количества нулей
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
    if (length(minimax[minimax != 1])>0)
    {
        minimax = minimax[minimax != 1]
    } else 
    {
        minimax = 1
    }
    return (minimax)
}

#функция-условие выбора слогов
convex = function (minmaxii,fh,d)
{
    convexi = matrix (0,(length(minmaxii)))
    
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

#функция-условие выбора слогов 2 
convexALT = function (minmaxii,fh,d)
{
    convexi = matrix (0,(length(minmaxii)))
    
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
        }
        #если предыдущая точка была мин и у нас макс (плюс проверка на два очень близких пика)
        if ((fh[minmaxi[i]]-fh[minmaxi[i-1]])>0.08 && flag == 0 )
        {
            steps = steps + 1
            convexi[steps] = minmaxi[i-1]
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
for (j in 1:1)
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
            #p = matrix(0,length(x)/2)
            if (st > 300 && st < 700 )
            {
                fh[i] = abs(x[k])
            }
            #fh[i] = max(p)
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
    b = convexALT(a,fh,0.2)
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
    #концы слов
    pause[t,1] = (words[2*j]*windowlen*x1@samp.rate)
    t = t+1
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

bb = convexALT(aa,fhii,0.2)
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



#------------------------------------------------------------
y4= Sys.time()
#10 MEL частот(12 опрных точек)
M = c(401.25, 622.50, 843.75, 1065.00, 1286.25, 1507.50, 1728.74, 1949.99, 2171.24, 2392.49, 2613.74, 2834.99)
#10 MEL частот в герцах(12 опрных точек)
h = c(300, 517.33, 781.90, 1103.97, 1496.04, 1973.32, 2554.33, 3261.62, 4122.63, 5170.76, 6446.70, 8000)


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

#на вход - слово, которое разбивается на фреймы wlen
MEL = function (signal,wlen)
{
    #нормируем
    signal = signal/max(abs(signal))
    len = length(signal)
    pos = 1
    #Число фреймов
    NumFrames = floor((len-wlen)/step) + 1
    MEL = matrix(0,NumFrames,11)
    
    #на спектре нашего фрейма
    f = floor((wlen+1)*h/x1@samp.rate)
    
    for (l in 1:NumFrames)
    {
        window = signal[pos:(pos+wlen-1)]
        pos = pos + 1
        S = matrix(0,10)
        for (m in 2:10)
        {
            for (k in 1:length(window))
            {
                S[m] = S[m]+ (abs(window[k])^2) * H(m,k,f)
            }
        }
        MEL[l,(2:11)] = S
        MEL[l,1] = sum (S^2)
    }
    return(MEL)
}


#MEL(signal[250000:350000],windowlen*x1@samp.rate)

pos = 1
wlen = 0.01*x1@samp.rate
step = wlen
len = length(signal)
#Разобьем все на фреймы, отсечем лишние частоты, полцчим melfh
NumFrames = floor((len-wlen)/step) + 1
melfh = matrix(0,length(signal))
for (i in 1:NumFrames)
{
    window = signal[pos:(pos+wlen-1)]
    pos = pos + step
    x = fft(window)
    bin = x1@samp.rate/length(x)
    st = 0
    p = matrix(0,length(x))
    for (k in 1:(length(x)/2))
    {
        st = st + bin
        if (st < 8000 )
        {
            p[k] = x[k]
        }
    }
    #print(wlen)
    #print (fft(p,inverse = T))
    p = Re(fft(p,inverse = T))
    melfh[((i-1)*wlen+1):((i)*wlen)] = p
}


# К каждому слогу применим mel-преоюразование
for (k in 2:nrow(pause))
{ 
    wsignal = melfh[floor(pause[k-1,1]):floor(pause[k,1])]
    if (sum(wsignal) != 0)
    {
        print(k)
    }
    if (sum(abs(wsignal))>epsilon)
    {
        wsignal = wsignal/max(wsignal)
        r=mel(wsignal)
    }
}

#----------------------
pos = 1
wlen = 0.01*x1@samp.rate
step = wlen
len = length(signal)
#Разобьем все на фреймы, отсечем лишние частоты, полцчим melfh
NumFrames = floor((len-wlen)/step) + 1
melfh = matrix(0,length(signal))
for (i in 1:NumFrames)
{
    window = signal[pos:(pos+wlen-1)]
    pos = pos + step
    x = fft(window)
    bin = x1@samp.rate/length(x)
    st = 0
    p = matrix(0,length(x))
    for (k in 1:(length(x)/2))
    {
        st = st + bin
        if (st < 8000 )
        {
            p[k] = x[k]
        }
    }
    #print(wlen)
    #print (fft(p,inverse = T))
    q = mel(p)
    #if (sum(abs(q))>0)
    #    print (q)
    p = Re(fft(p,inverse = T))
    melfh[((i-1)*wlen+1):((i)*wlen)] = p
}

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

#for (j in 1:(length(words)/2))
for (j in 1:1)
{
    #извлекаем нужные слоги
    #начало первого слога
    sbeg = nrow(pause)
    for (i in nrow(pause):1)
    {
        if (pause[i,1]>=(words[2*j-1]*windowlen*x1@samp.rate))
            sbeg = i
    }
    #Конец последнего
    send = nrow(pause)
    for (i in 1:nrow(pause))
    {
        if (pause[i,1]<=(words[2*j]*windowlen*x1@samp.rate) && pause[i,1] != pause[nrow(pause),1])
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
print (words)
