%Adaptive Modulation and Coding
%Baeline State

clc
%clear
close all
%% TDL/C Channel Parameters
SNR_c_dB=[-4.4 -1.2 -3.5 -5.2 -2.5 0 -2.2 -3.9 -7.4 -7.1 -10.7 -11.1 ...
	-5.1 -6.8 -8.7 -13.2 -13.9 -13.9 -15.8 -17.1 -16 -15.7 ...
	-21.6 -22.8];%power of maulipath components in TDL/C model

SNR_c=10.^(SNR_c_dB/10);

delay_c=[0 0.2099 .2219 .2329 .2176 .6366 .6448 .6560 .6584 .7935 ...
	.8213 .9336 1.2285 1.3083 2.1704 2.7105 4.2589 4.6003 5.4902 ...
	5.6077 6.3065 6.6374 7.0427 8.6523];%delay of multipath components in TDL/C model

rms_ds=20;%20ns rms delay spread of the channel 

time_unit=1000/40;

delay_cn=ceil(delay_c*rms_ds/time_unit);%normailized delay of multipath components

N_path=max(delay_cn+1);

del=0:max(delay_cn);

SNR_cn=zeros(size(del));

for i=1:length(del)
    for j=1:length(delay_cn)
        if delay_cn(j)==del(i)
            SNR_cn(i)=SNR_cn(i)+SNR_c(j);
        end
    end
end

SNR_cn=SNR_cn/max(SNR_cn);

SNR_cn=10*log10(SNR_cn);
%% Initial Parameters
M=[4 16 64];                %constellation size
    
N=1024;                     %number of subcarriers

N_pilot=N/8;                %total number of pilot samples

symb_d=1024/40;             %symbol duration

dp=del;                     %delay of multipath components in term of symbol duration

pl=SNR_cn;                  %logaritmic gain of multipath components  

N_path=length(dp)           %total number of paths

st=double(vpa(10/36,4));    %speed step size (equivalently 1km/h)

v=30*st;                    %initial velocity

v_ub=double(50*vpa(10/36,4));%upper band of velocity


v_lb=0;                     %lower band of velocity

ang=datasample(1:360,N_path);%intial vlaues of angle of arrivel in degree

f=3e9;                      %frequency

c=3e8;                      %seep of light

lambda=c/f;                 %wavelength

R=[1/2 1/2 2/3];                %code rate

C_s=[2 2 3];              %convelutional code states

%% other assessment parameters 
%SNR=datasample(6:0.25:26,1000);
%%
S=6:26;
SNR_counter=zeros(size(S));

N_be_Rayleigh=zeros(size(S));   %vector of nubmer of error blocks per SNR in Rayleigh Channel
N_be_AWGN=zeros(size(S));       %vector of number of error blocks per SNR in AWGN Channel

%%% SNR Levels %%%
%if 6<=SNR<12 ->QPSK
%elseif 12<=SNR<22 ->16QAM
%esleif 22<=SNR -> 64QAM

%SNR=datasample(6:0.25:26,1000);
SNR=10
%% AMC LOOP
tic
for k=1:length(SNR)
    if SNR(k)>=6
        SNR_round=floor(SNR(k));
        SNR_counter(SNR_round-5)=SNR_counter(SNR_round-5)+1;

        %% Modulation and Coding Choosing
        if SNR(k)>=6 && SNR(k)<12
            r=R(1);     %codaing rate
            m=M(1);      %modulation size
            c_s=C_s(1); %coding state
        elseif SNR(k)>=12 && SNR(k)<22
            r=R(2);
            m=M(2);
            c_s=C_s(2);
        elseif SNR(k)>=22
            r=R(3);
            m=M(3);
            c_s=C_s(3);
        end

        %% pilot generation
        decimal_pilot=datasample(0:m-1,N_pilot);        %decimal pilot samples
        pilot=qammod(decimal_pilot,m,'Gray');                                  %modulated pilot samples
           
        %% Bit Generation and Coding
        n_bit_message=N*7/8*log2(m)*r;      %number of message bits;
        b=randi([0 1],1,n_bit_message);     %message bit string
        b_coded=convcode(b,c_s);                %coded message bit string

        %%  OFDM Frame Generation
        s1=digqammod(b_coded, m);   %modulated bit strings
        s2=[];                      %OFDM Frame Before IFFT and Adding Cyclic Prefix
        for j=1:7:length(s1)-6
            s2=[s2 s1(j:j+6)  pilot(1+(j-1)/7)];
        end
        s2_ifft=ifft(s2);                   %OFDM Symbol After IFFT
        cp=s2_ifft(end-length(dp)+2:end);   %cylic prefix
        s=[cp s2_ifft];                     %Adding Cyclic Prefix to the Frame
    
        %% Rayleigh Channel
        %%%%%%%%%%%%%beginning of velocity%%%%%%%%%%%%%%%%
        if v==v_ub
            v=v+datasample([-1 0],1)*st;

        elseif v==v_lb
            v=v+datasample([0 1],1)*st;

        elseif v~=v_lb && v~=v_ub
            v=v+datasample([-1 0 1],1)*st;

        end

        %%%%%%%%%%%%%end of velocity%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%begging of angle of arrrivel%%%%%%%%
        ang_diff=0.5*datasample([-1 0 1],N_path);%difference of angle of arrivel in degree
        ang=ang+ang_diff;%new angel of arrive in degree
        ang=deg2rad(ang);%new angle of arrivel in degree

        %%%%%%%%%%%%%end of angle of arrivel%%%%%%%%%%%%%

        fd=v/lambda*cos(ang);%doppler frequency
        phi=rand(1,N_path)*2*pi;%random phase of arrived signal
        h=(cos(2*pi*fd+phi)+1j*sin(2*pi*fd+phi)).*(10.^(pl/10));
        if isa(h, 'sym')
            h = double(h);
        end
        
        %% Transmisson of Signal Through the Channel
        y=conv(s,h);                    %transmission of siganl through the rayleigh channel
        r=awgn(y,SNR(k),'measured');    %adding white noise with specific SNR

        %% Channel Estimation Using LS Algorithm
        r_N=r(length(cp)+1:length(cp)+length(s2_ifft)); %removing the cyclic prefix
        R_N=fft(r_N);                                   %fft of the received symbol block
        P_r=R_N(8:8:end);                               %speration of pilots of received block
        H_hat=P_r./pilot;                               %frequecny responce of the channel using LS Algorithm
        h1_hat=ifft(H_hat);                              %impulse responce of the channel
        h_hat=h1_hat(1:N_path);                          %considering the first 4 samples of the h_hat

        %% Equalization Using ZF Equalizer
        H_est=fft(h_hat,N);     %frequency Responce of Equlizer
        S_hat=R_N./(H_est);     %equlization
        S_hat(8:8:N)=[];        %removing the pilot symbols
        %% Detection for Rayleigh Channel
        s_hat=digqamdemod(S_hat,m);
        s_hat_decoded=convdecode(s_hat,c_s);
        err=biterr(s_hat_decoded,b);

        if err>0
            N_be_Rayleigh(SNR_round-5)=N_be_Rayleigh(SNR_round-5)+1;
        end

        %% Detection for AWGN Channel
        r2=awgn(s1,SNR(k),'measured');%trasmission through AWGN channel;
        s_hat2=digqamdemod(r2,m);
        s_hat2_decoded=convdecode(s_hat2,c_s);
        err2=biterr(s_hat2_decoded,b);
        
        if err2>0
            N_be_AWGN(SNR_round-5)=N_be_AWGN(SNR_round-5)+1;
        end
                
    end
end
toc

%% Throughput  and BLER Calculation
effect_blocks_R=SNR_counter-N_be_Rayleigh;
effect_blocks_A=SNR_counter-N_be_AWGN;

effect_bits_R=7/8*1024*[ones(1,length(6:11))*1/2*2 ones(1,length(12:21))*1/2*4 ...
                        ones(1,length(22:26))*2/3*6].*effect_blocks_R;

effect_bits_A=7/8*1024*[ones(1,length(6:11))*1/2*2 ones(1,length(12:21))*1/2*4 ...
                        ones(1,length(22:26))*2/3*6].*effect_blocks_A;

T_R=(effect_bits_R)./(symb_d*SNR_counter);

T_A=(effect_bits_A)./(symb_d*SNR_counter);

R_be_R=N_be_Rayleigh./SNR_counter;
R_be_A=N_be_AWGN./SNR_counter;

%% plot throughput
figure(1)
hold on

bar(S, T_R, 'BarWidth', 0.5, 'FaceAlpha', 0.6, 'EdgeColor', 'k','FaceColor','[0.5 0.5 .5]')
bar(S, T_A, 'BarWidth', 0.2, 'FaceAlpha', 0.6, 'EdgeColor', 'k','FaceColor','[ 0 0 0]')

ylabel('Throughput(Mbps)','Interpreter','latex','FontSize',14)
xlabel('$\gamma$(dB)','interpreter','latex','fontsize',14)

xticks(0:30);
xlim([4.5 26.5])
grid on

legend('Rayleigh','AWGN','interpreter','latex','fontsize',12,'location','northwest')
hold off
%% plot BLER
figure(2)
hold on

bar(S, R_be_R, 'BarWidth', 0.5, 'FaceAlpha',0.6 , 'EdgeColor', 'k','FaceColor','[0.5 0.5 0.5]')
bar(S, R_be_A, 'BarWidth', 0.2, 'FaceAlpha', 0.6, 'EdgeColor', 'k','FaceColor','[0 0 0]')

ylabel('BLER','Interpreter','latex','FontSize',14)
xlabel('$\gamma$(dB)','interpreter','latex','fontsize',14)

xticks(0:30);
xlim([4.5 26.5])
grid on

legend('Rayleigh','AWGN','interpreter','latex','fontsize',12,'location','northeast')


%% data saving
writematrix([R_be_A' R_be_R'],'D:\AMC_Baseline\BLER.xlsx')
writematrix([T_A' T_R'],'D:\AMC_Baseline\Throughput.xlsx')


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%# functions #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [r]=convdecode(x,k)

    if k==2 %1/2
    
        G=[56721 61713];
        %G=[133 171];
        %G=[4335 5723];
        %G=[561 753];
        %G=[447254 627324]
        %G=[716502 514576]
        L=15;
    
    elseif k==3 %R=2/3
    
        G=[63 15 46;
        32 65 61];
    
        L=[6 6];
    
    elseif k==4 %r=3/4
    
        G=[40 14 34 60;
        04 64 20 70;
        34 00 60 64];
    
        L=[6 6 6];
    end

    trellis=poly2trellis(L, G);

    %tb = ceil(length(x)/2);  
    tb=10*max(L);

    r=vitdec(x, trellis, tb, 'trunc', 'hard');

end

function [r]=convcode(x,k)

    if k==2
    
        G=[56721 61713];
        %G=[133 171];
        %G=[4335 5723];
        %G=[561 753];
        %G=[447254 627324]
        %G=[716502 514576]
        L=15;
    
    elseif k==3
    
        G=[63 15 46;
        32 65 61];
    
        L=[6 6];
    
    elseif k==4
    
        G=[40 14 34 60;
        04 64 20 70;
        34 00 60 64];
    
        L=[6 6 6];
    end
        trellis=poly2trellis(L, G);

    r=convenc(x, trellis);

end

function [y]=digqammod(x,N)
    y=[];
    if (N==4| N==16| N==64| N==256) & mod(length(x),log2(N))==0
        M=log2(N);
        for i=1:M:length(x)-M+1
            y=[y,bi2de(x(i:i+M-1),'left-msb')];
        end
        y=qammod(y,N,'Gray');
    end
end


function [x]=digqamdemod(y,N)
    if (N==4| N==16| N==64| N==256)
        x=qamdemod(y,N,'Gray');
        x=de2bi(x,log2(N),'left-msb');
        x=reshape(x',1,[]);
    end
end

