%% Training and dataset generation
function predictInitializations(modelPath,boxPath)
    % Generates predictions for the entire dataset and uses those for
    % initialization of unlabeled frames.

        if nargin < 1 || isempty(modelPath)
            modelPath = uibrowse([],[],'Select model folder...', 'dir');
            if isempty(modelPath) || ~exists(modelPath); return; end
        end

        % TODO: better system for choosing final vs best validation model
%         if exists(ff(modelPath, 'final_model.h5'))
%             numValidationSamples = numel(loadvar(ff(modelPath,'training_info.mat'),'val_idx'));
% %             numWeights = numel(dir_files(ff(modelPath,'weights')));
%             if numValidationSamples < 500
%                 modelPath = ff(modelPath,'final_model.h5');
%             end
%         end
        numValidationSamples = numel(loadvar(ff(modelPath,'training_info.mat'),'val_idx'));
        if exists(ff(modelPath, 'best_model.h5')) && numValidationSamples > 500
            modelPath = ff(modelPath, 'best_model.h5');
        else
            modelPath = ff(modelPath, 'final_model.h5');
        end

        % Predict
        preds = predict_box(boxPath, modelPath, false);

        % Save
        labels.initialization = preds.positions_pred;
        saveLabels();

        % Update status
        isInitialized = squeeze(all(~isnan(labels.initialization),2));
        numInitialized = sum(all(isInitialized,1));
        ui.status.framesInitialized.String = sprintf('Initialized: %d/%d (%.2f%%)', numInitialized, numFrames, numInitialized/numFrames*100);

        % Update status bars
        status = getStatus();
        ui.status.fullImg.CData = status;
        zoom_idx = ui.status.zoomImg.XData > 0 & ui.status.zoomImg.XData <= size(status,2);
        ui.status.zoomImg.CData(:,zoom_idx) = status(:,ui.status.zoomImg.XData(zoom_idx));

        % Log event
        addToHistory(['Initialized with model: ' modelPath])

        % Calculate error rate on labels
        labeled = all(getStatus() == 2,1);
        pos_gt = labels.positions(:,:,labeled);
        pos_pred = labels.initialization(:,:,labeled);
        pred_metrics = compute_errors(pos_pred,pos_gt);

        % Display errors
        printf('Error: mean = %.2f, s.d. = %.2f', mean(pred_metrics.euclidean(:)), std(pred_metrics.euclidean(:)))
        prcs = [50 75 90];
        prc_errs = prctile(pred_metrics.euclidean(:), prcs);
        for i = 1:numel(prcs)
            printf('       %d%% = %.3f', prcs(i), prc_errs(i))
        end
    end