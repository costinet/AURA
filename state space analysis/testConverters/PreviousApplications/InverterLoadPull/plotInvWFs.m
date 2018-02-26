function [ output_args ] = plotInvWFs( t, ys, R, colors, Ts)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    subplot(3,1,1);
    hold on;
    plot([t t+Ts/2], [ys(1,:), -ys(1,:)], 'Color', colors(1,:));
    plot([t t+Ts/2], [ys(3,:), -ys(3,:)], 'Color', colors(2,:));
    plot([t t+Ts/2], [ys(4,:)*R, -ys(4,:)*R], 'Color', colors(4,:));
    legend('V_p',  'v_{cr}',  'V_{out}');

    subplot(3,1,2);
    hold on;
    plot([t t+Ts/2], [ys(5,:), -ys(5,:)], 'Color', colors(2,:));
    legend(  'v_{cr}' );

    subplot(3,1,3);
    hold on;
    plot([t t+Ts/2], [ys(2,:), -ys(2,:)], 'Color', colors(1,:));
    plot([t t+Ts/2], [ys(4,:), -ys(4,:)], 'Color', colors(2,:));
    legend('i_{lr}',  'i_{lo}');

end

