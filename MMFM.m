function [state,R]=MMFM(pkts,a,N,M,npoints)

%Uso:
    %Uma vez modelada a fun��o de autocovari�ncia por uma exponencial do tipo
    %F(t)=var(pkts)*exp(-a*t), o programa utiliza o valor de 'a' encontrado e
    %calcula as intensidades de transi��es de estados correspondentes aos
    %nascimentos e mortes de uma cadeia de Markov do tipo "birth and death". 
    %O programa chama a fun��o 'birthdeath' e gera amostras segundo o MMFM (Modulated Markov Fluid Model)
%Par�metros:
    %pkts=S�rie de tr�fego a ser modelada
    %a=proveniente da modelagem da fun��o de autocovari�ncia
    %N=n�mero de fontes agregadas
    %npoints=tamanho da s�rie sint�tica gerada
    %M=A cadeia de markov possui M+1 estados. Variando de zero a M.
    %m=m�dia da s�rie
    %v=vari�ncia da s�rie
    %R=matriz de transi��es de estados
        %Inputs:
            %pkts,a,N,M,npoints
        %Outputs:
            %state, R        
    %Fl�vio Geraldo && C�lio Costa. UFG-Universidade Federal de Goi�s.
    %Trabalho realizado sob orienta��o do Prof. Dr. Fl�vio H. T. Vieira.
    
m=mean(pkts);
v=var(pkts);

%Modela as probabilidades do tr�fego atrav�s de uma distribui��o binomial
beta=a/(1+(N*((N*m)^2))/(M*N*v));
alfa=a-beta;
A=(N*v)/(N*m)+(N*m)/M;
%Inicializa a matriz de transi��es de estados
R=zeros(M+1);
%loop para cadeia de markov birth-death
M=M+1; %O programa considera estados {1,2,...,M+1}
for i=1:M
    if i<M
        R(i,i+1)=(M-i)*alfa;
    end
    if i>1
        R(i,i-1)=(i-1)*beta;
    end
    R(i,i)=0;
    for j=1:M
        if abs(i-j)>1
            R(i,j)=0;
        end
    end
end
%Nascimentos
nasc=M*alfa;
fl=[nasc:-alfa:alfa];
%Mortes
morte=M*beta;
fm=[beta:beta:morte];
[tjump, state] = birthdeath(npoints, fl, fm,A);
eixo=[1:npoints];
%plot(eixo,state);   %%%plota mudan�as de estados sem considerar tempos de
%transi��es
%figure;
[i,j]=hist(state,50);
y=i/npoints;
plot(j,y);
hold on
[p,k]=hist(pkts,50);
yy=p/npoints;
plot(k,yy,'r-');
media=mean(pkts);
tamanho=length(pkts);
poiss=poissrnd(media,tamanho,1);
[pp,kk]=hist(poiss,50);
yy=pp/npoints;
plot(kk,yy,'y-');
legend('MMFM','real','Poisson');
xlabel('Taxa de chegadas - pps');
ylabel('Probabilidade');
figure;
plot(tjump,state);
xlabel('tempo(segundos)');
ylabel('Taxa de chegadas - pps');

%Chama programa birthdeath
function [tjump, state] = birthdeath(npoints, flambda, fmu,A)
% BIRTHDEATH generate a trajectory of a birth-death process
%
% [tjump, state] = birthdeath(npoints[, flambda, fmu])
%
% Inputs: npoints - length of the trajectory
%         
% Outputs: tjump - jump times
%          state - states of the embedded Markov chain

 i=1;     %initial value, start on level i
 tjump(1)=0;  %start at time 0
 state(1)=i;  %at time 0: level i
 for k=2:npoints
    % compute the intensities
    lambda_i=flambda(i);
    mu_i=fmu(i);

    time=-log(rand)./(lambda_i+mu_i);      % Inter-step times:
                                       % Exp(lambda_i+mu_i)-distributed

    if rand<=lambda_i./(lambda_i+mu_i)
      i=i+1;     % birth
    else
      i=i-1;     % death
    end          %if
    state(k)=i;
    tjump(k)=time;
 end              %for i

 tjump=cumsum(tjump);     %cumulative jump times

 state=state-1;
 state=state*A;





