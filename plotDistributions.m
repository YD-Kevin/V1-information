function plotDistributions(r,n,pos)


    meanPoints=linspace(0,1,n+2);
    meanPoints=meanPoints(2:end-1);
    
    
    xRange=linspace(0,1,1000);
    figure('Color','white');
    hold on;
    xlim([0,1]);
    
    
    
    for i=1:n
       
        mu=meanPoints(i);
        pdfValues=normpdf(xRange,mu,r);
        
        
    
        plot(xRange,pdfValues,'LineWidth',1.2);
        
        
        posDensity=normpdf(pos,mu,r);
        plot([pos,pos],[0,posDensity],':','Color',[0.7,0.7,0.7],'LineWidth',0.8);
        plot(pos,posDensity,'ro','MarkerFaceColor','r','MarkerSize',8);
    end
    
    
    plot([pos,pos],[0,max(ylim)],'k:','LineWidth',3);
    
    
    
    title(sprintf('%d RFs (\\sigma=%.2f)',n,r));
    xlabel('x');
    ylabel('pdf');
    box on;
    grid on;


end