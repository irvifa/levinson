# Fungsi ini adalah penerapan dari algoritma Levinson
# Fungsi ini menerima 2 buah input
# r = matrix Toeplitz Simetris
# b = matrix yang memenuhi persamaan Tx = b
# Output berupa matrix x yang merupakan solusi

function x = levinson (r, b)
   # menyimpan nilai r0 dari matrix Toeplitz
   r0 = r(1);
   # Inisialisasi matrix r => [r1 r2 r3 ... rn-1]
   r  = (r(2:length(r))/r(1))';
   
   # Menghitung panjang matrix b 
   n  = length(b); 
   
   # Inisialisasi matrix v,x,y,z
   y  = zeros(n,1);
   v  = zeros(n,1); 
   z  = zeros(n,1); 
   x  = zeros(n,1);
   
   # Mengisi elemen tertentu pada matrix
   y(1) = -r(1);  
   x(1) = b(1); 
   
   # Inisialisasi nilai alfa dan beta
   beta = 1; 
   alfa = -r(1);
   
   # Looping
   for k=1:n-1
    # Inisialisasi nilai beta baru
    beta = (1-alfa^2)*beta;
    # Inisialisasi nilai mu
    mu   = (b(k+1)-r(1:k)'*x(k:-1:1))/beta;
    # Mengganti elemen pada matrix v dan x
    v(1:k)   = x(1:k)+ mu*y(k:-1:1);
    x(1:k+1) = [v(1:k);mu];
    if k<n-1
      # Inisialisasi nilai alfa baru
      alfa      = -(r(k+1)+r(1:k)'*y(k:-1:1))/beta;
      # Mengganti elemen pada matrix z
      z(1:k)    = y(1:k) + alfa*y(k:-1:1);
      y(1:k+1)  = [z(1:k);alfa];
    end
   end
   x = x/r0;
   x = x';
endfunction