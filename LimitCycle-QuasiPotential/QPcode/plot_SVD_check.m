function [] = plot_SVD_check(A)
        %A(1:nt,1:d,1:d);
        sa=size(A); nt = sa(1); d=min(sa(2:3));
        M=zeros(sa(2),sa(3));
        S = zeros(nt,d);
        for t = 1: nt 
            M(1:d,1:d)=A(t,1:d,1:d);
            S(t,1:d)= rcond(M); %svd(M);
        end
        for id = 1:d 
            plot([1:nt],log10(S(1:nt,id)));hold on ;grid on;
        end 
        ylabel('log10','FontSize',20)
        hold off
    end 