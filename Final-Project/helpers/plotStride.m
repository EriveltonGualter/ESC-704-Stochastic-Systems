function peaks = plotStride(varargin)
    A = varargin{1};
    if nargin == 2
        peaks = varargin{2};
    else
        a = detrend(A);
        [~,peaks] = findpeaks(a,'MinPeakProminence',(max(a)-min(a))/4);
    end
    %P         = mean(diff(peaks));

    for i=1:length(peaks)-1
        sa = A(peaks(i):1:peaks(i+1));
        plot(linspace(0,100,length(sa)), sa, '-k');
    end
    xlabel('% Stride');
    ylabel('Angle, rad');
end