function xdot = sir_eqn(x,t)
    % Parameter values
    beta=0.1;
    mu=0.05;

    % Define variables
    s = x(1);
    y = x(2);
    r = x(3);

    % Define ODEs
    ds=-beta*s*y;
    dy=beta*s*y-mu*y;
    dr=mu*y;

    % Return gradients
    xdot = [ds,dy,dr];
endfunction

t = linspace(0, 200, 2001)+.1;
x0=[0.99,0.01,0];

x = lsode("sir_eqn",x0, t);

out=[transpose(t),x];
save -ascii sir_octave.out out;

!head sir_octave.out

plot(t,x(:,1),"-r",t,x(:,2),"-g",t,x(:,3),"-b")
xlim([0 200])
xlabel("Time","fontweight","bold")
ylabel("Number","fontweight","bold")
h = legend("S","I","R");
legend(h,"show")
