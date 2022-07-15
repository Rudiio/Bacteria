clear all
close all

save=1

%% Donnée exp

bac=0; % 0 Coli, 1 pseudomonas

if bac==0
    divpara=0.046345606864877;
    nnew=2;
else
    divpara=0.072322790728674;
    nnew=3;
end

%Nom des fichiers
if bac==0 %Donnï¿½e toutes les 5 mins
    FILE={'coli/140930_MG1655_LB1-fld0';'coli/140930_MG1655_LB1-fld1';'coli/140930_MG1655_LB1-fld3';'coli/141001_MG1655_LB1-fld0';'coli/141001_MG1655_LB1-fld2';'coli/141001_MG1655_LB1-fld3';'coli/141003_MG1655_LB1-fld1'};
elseif bac==1 %Donnï¿½e toutes les 10 mins
    FILE={'pseudomonas/100806_pao1';'pseudomonas/101105_pao1_bottom';'pseudomonas/101105_pao1_top';'pseudomonas/101108_pao1';'pseudomonas/101110_pao1';'pseudomonas/101112_pao1';'pseudomonas/101114_pao1';'pseudomonas/101115_pao1';'pseudomonas/101116_pao1';'pseudomonas/101119_pao1'};
end

%Données que l'on veut sauvegarder
NumberBacE=zeros(length(FILE),100); %Nombre de bactérie au cours du temps
IncrE=[];
LE=[];
LbirthE=[];
LdivE=[];
GRE=[];
GenE=[]; % Sauvegarde des générations
Firstdiv=[];
ThetaE=[];
lifespanE=[]

%Extraction des données
for cas=1:length(FILE)
    file_path=FILE{cas};
   
    [cname, frame, time, ncells, schnitz, birth, lifespan, cell_length, cell_area, x_cell, y_cell, angle, ColoCenterDist, ColoBorderDist, AverageColoBorderDist, grate, grcoef, grateinst, width,...
    totpix_fluo, siderosd, meanpix_fluo,siderobg,mChetot,mChesd,mChemed,mChebg,siderongb,areangb,siderongbCell,areangbCell,siderongbl,areangbl,siderongblCell,areangblCell,AvgNgbAngle, correltot,correlsd,...
    correlmed, Old_pole_X, Old_pole_Y, Young_pole_X, Young_pole_Y, OldPoleAge, Bary_Neigh_X, Bary_Neigh_Y]=textread([file_path '.txt'],'%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f'...
    ,-1,'headerlines',1,'delimiter','\t','endofline','\r\n');

%     for k=9:max(birth)
%         Ind=find(birth==k); %Toutes ces cellules sont nées aux memes moments
%         for i1=1:length(Ind)
%             for i2=1:length(Ind)
%                 if (i1<i2)&&(strcmp(cname{Ind(i1)}(1:end-1),cname{Ind(i2)}(1:end-1)))
%                     ind1=find(schnitz==schnitz(Ind(i1)));
%                     ind2=find(schnitz==schnitz(Ind(i2)));
%                     ThetaE=[ThetaE angle(ind1(1))-angle(ind2(1))];
%                 end
%             end
%         end
%     end

   for fr=2:max(frame)
       NumberBacE(cas,fr)=length(find(frame==fr));
   end

    for i=2:max(schnitz)
        Ind=find(schnitz==i);
        if (frame(Ind(end))~=max(frame))
            lifespanE= [lifespanE lifespan(Ind(1))];
        end
        if (length(Ind)>1)&&(length(cname{Ind(1),1})>1)&&(frame(Ind(end))~=max(frame)) 
            %on exclut les cas ou le vecteur Ind est vide ou de longueur 1 (pbl dans la segmentation je pense)
            %On exclut le cas "artificiel" du debut (cname='0') et les cas 'OT' et 'OH' car les bacteries ne venait pas forcement de naitre au debut
            %on exclut le cas ou la bacterie est prï¿½sent au dernier temps i.e. elle ne s'est pas divisï¿½ avant la fin de l'experience
            IncrE=[IncrE cell_length(Ind(end))-cell_length(Ind(1))]; % calcul de l'increment (taux ï¿½ la division - taux a la naissance)
            LE=[LE cell_length(Ind)'];
            LdivE=[LdivE cell_length(Ind(end))];
            LbirthE=[LbirthE cell_length(Ind(1))]; %ajout Marie de la taille Ã  la naissance
            GRE=[GRE log(cell_length(Ind(end))/cell_length(Ind(1)))/(double(time(Ind(end)))-double(time(Ind(1))))];
            GenE=[GenE length(cname{Ind(1),1})];
            
        end
        if (length(cname{Ind(1),1})==3)
          Firstdiv=[Firstdiv cell_length(Ind(1))];
        end
    end
end
ThetaE=ThetaE-round(ThetaE/pi)*pi;

%% Données simulations

CAS=[1];
Ncas=length(CAS);
Nini=5;
r0=0.5;

%Données que l'on veut sauvegarder
NumberBacS=zeros(Ncas,Nini,300); %Nombre de bactérie au cours du temps
IncrS=[];
LS=[];
LbirthS=[];
LdivS=[];
GRS=[];
ThetaS=[];

for icas=1:Ncas
    cas=CAS(icas);
    
    for ini=1:Nini
        file=[ 'new' num2str(nnew ) '/cas=' num2str(cas) '_ini=' num2str(ini) '.dat'];
        densit1=fopen(file);
        A=fscanf(densit1,'%g %g %g %g',[8 inf]);
        Nb=A(1,:);
        Time=A(2,:);
        XX=A(3,:);
        YY=A(4,:);
        T=A(5,:);
        L1=A(6,:);
        L2=A(7,:);
        Sister=A(8,:);
        fclose(densit1);

        Tmax=Time(end);
        Nmax=length(T);
        NN=max(Nb);
        X=[XX' YY'];
        T=T-2*floor(T/(2*pi))*pi;
        L=L1+L2;
        
%         for k=1:max(Sister)
%             Ind=find(Sister==k);
%             if length(Ind)==2
%                 ind1=find((Nb==Nb(Ind(1))).*(mod(Time,3)==0));
%                 ind2=find((Nb==Nb(Ind(2))).*(mod(Time,3)==0));
%                 ThetaS=[ThetaS T(ind1(1))-T(ind2(1))];
%             end
%         end

        for it=1:Tmax
            NumberBacS(icas,ini,it)=length(find(Time==it));
        end

        for i=1:NN
            Ind=find(Nb==i);
            if length(Ind>0)
                if (length(Ind)>1)&&(Time(Ind(1))<Time(Ind(end)))&&(Time(Ind(1))>1)&&(Time(Ind(end))<Tmax)
                    IncrS=[IncrS L(Ind(end))-L(Ind(1))];
                    LS=[LS L(Ind)+2*r0];
                    LbirthS=[LbirthS L(Ind(1))+2*r0];
                    LdivS=[LdivS L(Ind(end))+2*r0];
                    GRS=[GRS 1/(Time(Ind(end))-Time(Ind(1)))*log(L(Ind(end))/L(Ind(1)))];
                end
            end
        end
    end
end
ThetaS=ThetaS-round(ThetaS/pi)*pi;


Incr_corr=IncrE(find(IncrE>0)); %on enlÃ¨ve les indices oÃ¹ l'incrÃ©ment est nÃ©gatif
Ldiv_corr=LdivE(find(IncrE>0));
[f_incr x_incr h]=ksdensity(Incr_corr,'Weights',Ldiv_corr,'support',[0 max(Incr_corr)*1.5]);

%% Plots



figure(1)
clf
subplot(1,2,1)
hold on
for cas=1:length(FILE)
    plot((2:3:(3*length(NumberBacE(cas,:))))-2,NumberBacE(cas,:),'-r')
end
for cas=1:Ncas
    for ini=1:Nini
        plot(1:300,squeeze(NumberBacS(cas,ini,:)),'-b')
    end
end
subplot(1,2,2)
hold on
for cas=1:length(FILE)
    plot((2:3:(3*length(NumberBacE(cas,:))))-2,NumberBacE(cas,:),'-r')
end
for cas=1:Ncas
    for ini=1:Nini
        plot((60:300)-59,squeeze(NumberBacS(cas,ini,60:end)),'-b')
    end
end
title('Number of bacteria as a function of time')
saveas(gcf,['Nb.png'])

figure(2)
clf
hold on
[f,xi,bw] =ksdensity(IncrE);
plot(xi,f,'-r','LineWidth',3)
[f,xi,bw] =ksdensity(IncrS/divpara);
plot(xi,f,'-b','LineWidth',3)
title('Distibution Incr')
legend(['Experimental data'],[' Simulations data'])
saveas(gcf,['Incr_s=' num2str(save) '.png'])
    
figure(3)
clf
hold on
[f,xi,bw] =ksdensity(LE);
plot(xi,f,'-r','LineWidth',3)
[f,xi,bw] =ksdensity(LS/divpara);
plot(xi,f,'-b','LineWidth',3)
title('Distibution length')
legend(['Experimental data'],[' Simulations data'])
saveas(gcf,['Length_s=' num2str(save) '.png'])

figure(4)
clf
hold on
[f,xi,bw] =ksdensity(LbirthE);
plot(xi,f,'-r','LineWidth',3)
[f,xi,bw] =ksdensity(LbirthS/divpara);
plot(xi,f,'-b','LineWidth',3)
title('Distibution birth length')
legend(['Experimental data'],[' Simulations data'])
saveas(gcf,['LengthB_s=' num2str(save) '.png'])

figure(5)
clf
hold on
[f,xi,bw] =ksdensity(LdivE);
plot(xi,f,'-r','LineWidth',3)
[f,xi,bw] =ksdensity(LdivS/divpara);
plot(xi,f,'-b','LineWidth',3)
title('Distibution division length')
legend(['Experimental data'],[' Simulations data'])
saveas(gcf,['LengthD_s=' num2str(save) '.png'])

GRE=GRE(find(GRE>=0));
figure(6)
clf
hold on
[f,xi,bw] =ksdensity(GRE);
plot(xi,f,'-r','LineWidth',3)
[f,xi,bw] =ksdensity(GRS);
plot(xi,f,'-b','LineWidth',3)
title('Distibution Growth rate')
legend(['Experimental data'],[' Simulations data'])
saveas(gcf,['GR_s=' num2str(save) '.png'])

figure(7)
clf
hold on
[f,xi,bw] =ksdensity(ThetaE);
plot(xi,f,'-r','LineWidth',3)
[f,xi,bw] =ksdensity(ThetaS);
plot(xi,f,'-b','LineWidth',3)
title('Distibution Angle at division')
legend(['Experimental data'],[' Simulations data'])
saveas(gcf,['dT_s=' num2str(save) '.png'])






