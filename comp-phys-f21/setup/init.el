;;  Set initial position of Emacs window ;
(setq initial-frame-alist '((top . 10) (left . 30)
			    (width . 100) (height . 40)))

;;  Frame appearence
(setq default-frame-alist
      '((top . 20) (left . 40)
	(width . 120) (height . 40)
	(cursor-color .     "blue")
	(foreground-color . "black")
	(background-color . "white")))

;; 16-Aug-2018: need to set these twice on Carbonate - have to use hex here
(add-to-list 'default-frame-alist '(foreground-color . "#000000"))
(add-to-list 'default-frame-alist '(background-color . "#ffffff"))

;; remove split screen on startup (08-Aug-2018)
(setq inhibit-startup-screen t)

;; set the format for the title bars (and thus the icon text) for 
;; each window			
(setq frame-title-format "Emacs - %f")
(setq icon-title-format "Emacs - %b")

;;  Key Bindings			
(global-set-key [C-f4]  'kill-this-buffer)
(global-set-key [M-f5]  'goto-line)

;; Text mode				
(add-hook 'text-mode 'auto-fill-mode)
(setq-default fill-column 100)

;; Default font (08-Aug-2018)
(set-default-font "DejaVu Sans Mono-14")

;; Fontification
(cond ((fboundp 'global-font-lock-mode)
            ;; Customize face attributes
            (setq font-lock-face-attributes
                  ;; Symbol-for-Face Foreground Background Bold Italic Underline
                  '((font-lock-comment-face       "DarkGreen")
                    (font-lock-string-face        "Sienna")
                    (font-lock-keyword-face       "Red")
                    (font-lock-function-name-face "Blue")
                    (font-lock-variable-name-face "SlateBlue")
                    (font-lock-type-face          "Maroon")
;;                    (font-lock-reference-face     "Purple")
                    ))
            ;; Load the font-lock package.
            (require 'font-lock)
            ;; Maximum colors
            (setq font-lock-maximum-decoration t)
            ;; Turn on font-lock in all modes that support it
            (global-font-lock-mode t)))			

;; (setq font-lock-support-mode 'lazy-lock-mode) ; 18-Oct-07 remove for java
(show-paren-mode 1)                              ; parantheses matching
(transient-mark-mode t)                          ; highlight mark region

;; Calendar/Time support
(setq european-calendar-style 't)
(setq calendar-week-start-day 1)
(display-time)

;; C support				
(load "cc-mode")
(add-hook 'c-mode 'auto-fill-mode)
(add-hook 'c++-mode 'auto-fill-mode)
(add-hook 'java-mode 'auto-fill-mode)
(setq auto-mode-alist
      (append
       '(("\\.C$"    . c++-mode)
	 ("\\.H$"    . c++-mode)
	 ("\\.cc$"   . c++-mode)
	 ("\\.cpp$"  . c++-mode)
	 ("\\.hh$"   . c++-mode)
	 ("\\.hpp$"  . c++-mode)
	 ("\\.c$"    . c-mode)
	 ("\\.sc$"   . c-mode)
	 ("\\.h$"    . c-mode)
	 ("\\.m$"    . objc-mode)
	 ("\\.java$" . java-mode)
	 ("\\.sa$"   . asm-mode)
	 ) auto-mode-alist))
