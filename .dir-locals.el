;; Local variables specific to this project.
;;
;; The add-hook functionality requires that the following code is included in the main dotfile:
;; ;; Provide a new MAJORMODE-local-vars-hook
;; ;; reference: https://emacs.stackexchange.com/a/29662
;; (add-hook 'hack-local-variables-hook 'run-local-vars-mode-hook)
;; (defun run-local-vars-mode-hook ()
;;   "Run a hook for the major-mode after the local variables have been processed."
;;   (run-hooks (intern (concat (symbol-name major-mode) "-local-vars-hook"))))
;;
;; usage:
;; ((nil . ((eval . (add-hook 'c++-mode-local-vars-hook  (lambda () (c-set-offset 'innamespace '+))))
;;          )
;;       )
;;  )
;;
;; Reference:
;; https://www.gnu.org/software/emacs/manual/html_node/emacs/Directory-Variables.html

((nil . ((indent-tabs-mode . t)
         (c-basic-offset . 8)
         (tab-width . 8)
         (fill-column . 120)
         (eval . (add-hook 'c-mode-local-vars-hook  (lambda () (c-set-offset 'topmost-intro-cont '+))))
         (eval . (add-hook 'c-mode-local-vars-hook  (lambda () (modify-syntax-entry ?_ "w")))) ; Make _ symbols count as part of words in c-mode.
         )
      )
 (c-mode (helm-make-build-dir . "build_debug_2D"))
 )
