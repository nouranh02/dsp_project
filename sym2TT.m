function x = sym2TT(S)
    % Define the DTMF frequency pairs
    F = [697, 770, 852, 941; 1209, 1336, 1477, 1633];
    
    % Define the DTMF symbol table
    DTMF = ['1', '2', '3', 'A'; '4', '5', '6', 'B'; '7', '8', '9', 'C'; '*', '0', '#', 'D'];
    
    % Find the row and column indices of the input symbol in the DTMF table
    [r, c] = find(DTMF == S);
    
    % Generate the corresponding DTMF signal
    f1 = F(1, r);
    f2 = F(2, c);
    fs = 8000;
    Ts = 1/fs;
    t = 0:Ts:0.1;
    x = sin(2*pi*f1*t) + sin(2*pi*f2*t);
end