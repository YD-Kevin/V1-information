function sd = pdfStats(x, pdf_values,pos)


x = x(:);                
pdf_values = pdf_values(:); 


norm_factor = trapz(x, pdf_values);
if norm_factor <= 0
    error('概率密度积分非正，无效输入');
end
pdf_normalized = pdf_values / norm_factor;


mu = trapz(x, x .* pdf_normalized);

Ex2 = trapz(x, (x.^2) .* pdf_normalized);


sigma2 = Ex2 -2*mu*pos + pos^2;


if sigma2 < 0
    sigma2 = abs(sigma2); 
end

sd=sqrt(sigma2);
end