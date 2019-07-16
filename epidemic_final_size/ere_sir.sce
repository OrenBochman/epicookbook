function final_size=SIR_fs(N,bet,gamm)
    final_size = zeros(N+1,1);
    final_size(2) = 1;
    for Z2 = 0:N;
        for Z1 = Z2+1:N-1
            p1 = 1 / ( 1 + gamm/(bet*(N-Z1)));
            final_size(Z1+2) = final_size(Z1+2) + final_size(Z1+1)*p1;       
            final_size(Z1+1) = final_size(Z1+1)*(1-p1);
        end
    end
endfunction

N = 20;                       
bet = 2/(N-1);
gamm = 1;

tic();
final_size=SIR_fs(N,bet,gamm);
toc()

bar(0:N,final_size)
