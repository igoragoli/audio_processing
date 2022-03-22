function [alphaloqo,yloqo,bbb] = uncontretousoutil(donnees,nn1,nn2,classes,sig2,nu, plot);

use_octave=1;

ttloqo=clock;

cl1_d1=donnees(1,1:nn1);
cl1_d2=donnees(2,1:nn1);
cl2_d1=donnees(1,nn1+1:end);
cl2_d2=donnees(2,nn1+1:end);

mind1 = min([cl1_d1 cl2_d1]);
maxd1 = max([cl1_d1 cl2_d1]);
mind2 = min([cl1_d2 cl2_d2]);
maxd2 = max([cl1_d2 cl2_d2]);

for ii=1:nn1+nn2
  if (ii<=nn1)
    vaii1=cl1_d1(ii);
    vaii2=cl1_d2(ii);
  else
    vaii1=cl2_d1(ii-nn1);
    vaii2=cl2_d2(ii-nn1);
  end;
  for jj=1:nn1+nn2
    if (jj<=nn1)
      vajj1=cl1_d1(jj);
      vajj2=cl1_d2(jj);
    else
      vajj1=cl2_d1(jj-nn1);
      vajj2=cl2_d2(jj-nn1);
    end;

    % noyau RBF gaussien
    matsca(ii,jj) = lagis_rbf_gaussien([vaii1 vaii2]', [vajj1 vajj2]', sig2);
  end;
end;
ppx=-ones(1,nn1);
ppy= ones(1,nn2);
matpp=[ppx ppy]'*[ppx ppy];

c = -ones(nn1+nn2,1);
H = matpp.*matsca;
A = [ppx ppy];
b = 0;
l = zeros(nn1+nn2,1);
u = 1/(nu*(nn1+nn2))*ones(nn1+nn2,1);
%u = 16.0*ones(2*nn,1);   %%% 'nu' est contrôlé ici
[alphaloqo,yloqo] = pr_loqo3(c, H, A, b, l, u, use_octave);
eeloqo=etime(clock,ttloqo);

fprintf(1,'temps requis pour LOQO : %f\n',eeloqo);

alphaxloqo = alphaloqo(1:nn1)';        % alphas pour la première classe
alphayloqo = alphaloqo(nn1+1:end)';   % alphas pour la deuxième classe

%%%%%%%%%%%%%%%%%%%%%
alpha=alphaloqo;
alphax=alphaxloqo;
alphay=alphayloqo;

%%% calcul de la solution : calcul de w

ww1r = sum(-alphax.*cl1_d1) + sum(alphay.*cl2_d1);
ww2r = sum(-alphax.*cl1_d2) + sum(alphay.*cl2_d2);
normal = mean([ww1r ww2r]);
ww1=ww1r/normal;      % première dimension de w
ww2=ww2r/normal;      % deuxième dimension de w

bbbopt=1/(sqrt(ww1r^2 + ww2r^2));

tmpv2 = alpha.*classes';
for ii=1:nn1
  cl1_v(ii) = sign( tmpv2'*lagis_rbf_gaussien([cl1_d1(ii) cl1_d2(ii)]', donnees, sig2)' + bbbopt );
end;
for ii=1:nn2
  cl2_v(ii) = sign( tmpv2'*lagis_rbf_gaussien([cl2_d1(ii) cl2_d2(ii)]', donnees, sig2)' + bbbopt );
end;
bbb=bbbopt;

if plot
    figure();
    clf;
    subplot(211);
    title('sortie du SVM pour les points de la classe 1');
    grid on;
    hold on;
    plot(cl1_v,'r');
    hold off;
    subplot(212);
    title('sortie du SVM pour les points de la classe 2');
    grid on;
    hold on;
    plot(cl2_v,'b');
    hold off;
end


%%% au cas où
if bbb==-10000000
  fprintf(1,'b n a pas pu etre calcule correctement\n');
  bbb=-1.0
end;


%%% plot pour la route - LOQO
calcul_loqo=1;
if (calcul_loqo==1) && plot
  vsloqocl1 = alphaxloqo>max(alphaloqo)/20.0;
  vsloqocl2 = alphayloqo>max(alphaloqo)/20.0;
  figure();
  clf;
  grid on;
  hold on;
  title('les points magentas sont les vecteurs supports trouves par LOQO')
  plot(cl1_d1,cl1_d2,'*r');
  plot(cl2_d1,cl2_d2,'ob');
  for ii=1:nn1
    if (vsloqocl1(ii)>0)
      plot(cl1_d1(ii),cl1_d2(ii),'+m');
    end
  end
  for ii=1:nn2
    if (vsloqocl2(ii)>0)
      plot(cl2_d1(ii),cl2_d2(ii),'+m');
    end
  end
  xlabel('dimension 1');
  ylabel('dimension 2');
  hold off;
end


%%% calcul et plot de la fonction de décision obtenue
nzz=500;

mind1c=mind1;
maxd1c=maxd1;
pasd1=(maxd1c-mind1c)/nzz;

mind2c=mind2;
maxd2c=maxd2;
pasd2=(maxd2c-mind2c)/nzz;

ppd1res = [];
ppd2res = [];
ppd1=mind1c;
for ii=1:nzz
  ppd2=mind2c;
  for jj=1:nzz
     clnew(ii,jj) = sign( tmpv2'*lagis_rbf_gaussien([ppd1 ppd2]', donnees, sig2)' + bbbopt );
     if ii==1
       ppd2res=[ppd2res ppd2];
     end;
    ppd2=ppd2+pasd2;
  end;
  ppd1res=[ppd1res ppd1];
  ppd1=ppd1+pasd1;
end;

%%% plot pour la route
if plot
    figure();
    clf;
    title('frontiere entre les 2 classes');
    hold on;
    imagesc(ppd1res,ppd2res,-clnew');
    % attention : "mesh" inverse l'axe des abscisses et l'axe des ordonnées
    % => même comportement qu'avec "matlab", donc :-(
    %mesh(ppd2res,ppd1res,clnew);
    grid on;
    hold on;
    plot(cl1_d1,cl1_d2,'*r');
    plot(cl2_d1,cl2_d2,'ob');
    xlabel('dimension 1');
    ylabel('dimension 2');
    xlim([min(ppd1res) max(ppd1res)]);
    ylim([min(ppd2res) max(ppd2res)]);
    hold off;
    colormap('cool');
end

