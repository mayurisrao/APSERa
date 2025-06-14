 function diffsum = chromaticity1(fit_beam_der)
    % fit_beam_der=Gz;
    n_theta = size(fit_beam_der,1);
    N_phi = size(fit_beam_der,2);
    N_freq = size(fit_beam_der,3);

    beam_der_fit_1 = zeros(n_theta, N_phi, N_freq);
    beam_der_fit_1=fit_beam_der;
    lsum = 0;
    lsum1 = 0;
    diffsum = 0;

    for i = 1:(n_theta)
        for j = 1:N_phi
            l = squeeze(beam_der_fit_1(i, j, :)); % Extract 1D array
            f = find(diff(sign(l))); % Find zero crossings
          
            if ~isempty(f)
                temp_split = mat2cell(l, diff([0; f(:); numel(l)])); % Split the array
                for d = 1:length(temp_split)
                    segment = abs(temp_split{d});
                    lmax = max(segment);
                    lmin = min(segment);
                    lsum1 = lsum1 + lmax + lmin;
                    diff_val = lmax - lmin;
                    diffsum = diffsum + diff_val;
                    lsum = lsum + lmax;
                end
            else
                lmax = max(abs(l));
                lmin = min(abs(l));
                lsum1 = lsum1 + lmax + lmin;
                lsum = lsum + lmax;
                diffsum = diffsum + (lmax - lmin);
            end
        end
    end
 end
