clear all; close all; clc

Af = .9335; Bf = .0918; 
As = .9978;  Bs = .0143 
Bv(1:2) = [.48 .1]; 

nBL = 200; nPert = 400; nWO = 400;
nStrides = nBL + nPert + nWO;

%let's make error a visuomotor rotation (or clamp)
split = [zeros(nBL,1); ones(nPert,1); zeros(nWO,1)];
vis = [zeros(nBL,1); ones(nPert,1); zeros(nWO,1)];

[xf(1), xs(1), xv(1), x(1), e(1)] = deal(0);

for i=1:2
       
    for n=1:nStrides-1
        
        if split(n)==1
            pert(n) = 1;
        else
            pert(n) = 0;
        end
        e(n) = pert(n) - x(n);
        
        xf(n+1) = Af*xf(n) + Bf*e(n);
        
        xs(n+1) = As*xs(n) + Bs*e(n);
        
        x(n+1) = xf(n+1) + xs(n+1);
        
        
        
        if vis(n) == 1
            if i==1
                xv(n+1) = Bv(1)*e(n);
            elseif i==2
                xv(n+1) = Bv(2)*e(n);
            end
            xmotor(n) = x(n) + xv(n);
        else
            xv(n+1) = 0;
            xmotor(n) = x(n);
        end
        
    end
    
    subplot(1,2,i); hold on
    set(gcf, 'units','inches','pos',[5 5 17.5 5])
    xlabel('Stride number','fontsize',16)
    ylabel('Symmetry change','fontsize',16)
    ylim([-.5 1.2])
    
    if i==1
        plot(xmotor,'b','linewidth',4)
        title('Healthy','fontsize',18)
    else
        plot(xmotor,'color',[1 .647 0],'linewidth',4)
        title('Stroke survivor','fontsize',18)
    end
    plot(x,'r','linewidth',3)
    plot(xs,'g','linewidth',2)
    plot(xf,'c','linewidth',2)
    plot(xv,'m','linewidth',2)
    plot(split,'--','linewidth',2,'color',[0.5 0.5 0.5])
    legend({'motor output','adaptation','slow process','fast process','voluntary correction','perturbation'},'fontsize',12)
    legend boxoff
end
%print(gcf,'simulations','-r300','-dtiff')







