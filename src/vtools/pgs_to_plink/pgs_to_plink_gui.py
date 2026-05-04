"""
PGS → PLINK Converter GUI Wrapper
Оборачивает pgs_to_plink.py в простой tkinter интерфейс.
Все input() вызовы перехватываются через очередь.
"""

import tkinter as tk
from tkinter import ttk, scrolledtext, filedialog, messagebox
import threading
import queue
import sys
import os
import io

# ─── Путь к основному скрипту ────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(SCRIPT_DIR)

# ─── Перехват stdin / stdout / stderr ────────────────────────────────────────

class GuiOutput(io.TextIOBase):
    """Перенаправляет print() в виджет ScrolledText."""
    def __init__(self, widget: scrolledtext.ScrolledText):
        self.widget = widget

    def write(self, text: str) -> int:
        if text:
            self.widget.after(0, self._insert, text)
        return len(text)

    def _insert(self, text: str):
        self.widget.configure(state="normal")
        self.widget.insert(tk.END, text)
        self.widget.see(tk.END)
        self.widget.configure(state="disabled")

    def flush(self):
        pass


class GuiInput(io.TextIOBase):
    """
    Перехватывает input().
    Когда скрипт вызывает input(prompt), GUI показывает prompt в консоли
    и разблокирует поле ввода. Пользователь вводит ответ и жмёт Enter.
    Результат возвращается скрипту.
    """
    def __init__(self, output_widget: scrolledtext.ScrolledText,
                 entry_widget: tk.Entry, send_btn: tk.Button,
                 input_queue: queue.Queue):
        self.out = output_widget
        self.entry = entry_widget
        self.btn = send_btn
        self.q = input_queue

    def readline(self) -> str:
        # разблокируем поле ввода в главном потоке
        self.entry.after(0, self._unlock)
        value = self.q.get()          # блокируемся до ввода
        self.entry.after(0, self._lock)
        # эхо введённого значения в консоль
        self.out.after(0, self._echo, value)
        return value + "\n"

    def _unlock(self):
        self.entry.configure(state="normal")
        self.btn.configure(state="normal")
        self.entry.focus_set()

    def _lock(self):
        self.entry.configure(state="disabled")
        self.btn.configure(state="disabled")

    def _echo(self, value: str):
        self.out.configure(state="normal")
        self.out.insert(tk.END, value + "\n")
        self.out.see(tk.END)
        self.out.configure(state="disabled")

    def flush(self):
        pass


# ─── Основное окно ────────────────────────────────────────────────────────────

class PgsGuiApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("🧬 PGS → PLINK Converter v5.0")
        self.geometry("900x650")
        self.minsize(700, 500)
        self.configure(bg="#1e1e2e")

        self._input_queue: queue.Queue = queue.Queue()
        self._script_thread: threading.Thread | None = None

        self._build_ui()
        self._setup_io()

        # запускаем скрипт сразу при открытии
        self.after(200, self._run_script)

    # ── UI ──────────────────────────────────────────────────────────────────

    def _build_ui(self):
        # ── верхняя панель (заголовок + кнопки) ─────────────────────────────
        top = tk.Frame(self, bg="#181825", pady=6, padx=10)
        top.pack(fill="x")

        tk.Label(
            top, text="🧬  PGS → PLINK Converter", font=("Segoe UI", 14, "bold"),
            bg="#181825", fg="#cdd6f4"
        ).pack(side="left")

        self.restart_btn = tk.Button(
            top, text="⟳  Restart", command=self._restart,
            bg="#313244", fg="#cdd6f4", activebackground="#45475a",
            relief="flat", padx=10, pady=4, cursor="hand2",
            font=("Segoe UI", 10)
        )
        self.restart_btn.pack(side="right", padx=4)

        self.browse_btn = tk.Button(
            top, text="📂  Browse file…", command=self._browse_file,
            bg="#313244", fg="#cdd6f4", activebackground="#45475a",
            relief="flat", padx=10, pady=4, cursor="hand2",
            font=("Segoe UI", 10)
        )
        self.browse_btn.pack(side="right", padx=4)

        # ── консоль вывода ────────────────────────────────────────────────
        console_frame = tk.Frame(self, bg="#1e1e2e", padx=10, pady=(6, 0))
        console_frame.pack(fill="both", expand=True)

        self.console = scrolledtext.ScrolledText(
            console_frame,
            bg="#11111b", fg="#cdd6f4",
            insertbackground="#f5c2e7",
            font=("Consolas", 10),
            state="disabled",
            relief="flat",
            wrap="word",
            padx=8, pady=6
        )
        self.console.pack(fill="both", expand=True)

        # цветовые теги
        self.console.tag_config("warn",  foreground="#fab387")
        self.console.tag_config("error", foreground="#f38ba8")
        self.console.tag_config("ok",    foreground="#a6e3a1")
        self.console.tag_config("info",  foreground="#89dceb")

        # ── строка ввода ──────────────────────────────────────────────────
        input_frame = tk.Frame(self, bg="#1e1e2e", padx=10, pady=8)
        input_frame.pack(fill="x")

        prompt_lbl = tk.Label(
            input_frame, text="›", bg="#1e1e2e", fg="#89b4fa",
            font=("Consolas", 14, "bold")
        )
        prompt_lbl.pack(side="left", padx=(0, 4))

        self.entry = tk.Entry(
            input_frame,
            bg="#313244", fg="#cdd6f4",
            insertbackground="#f5c2e7",
            font=("Consolas", 11),
            relief="flat",
            state="disabled"
        )
        self.entry.pack(side="left", fill="x", expand=True, ipady=5)
        self.entry.bind("<Return>", self._send_input)

        self.send_btn = tk.Button(
            input_frame, text="Send", command=self._send_input,
            bg="#89b4fa", fg="#1e1e2e", activebackground="#b4befe",
            relief="flat", padx=14, pady=5, cursor="hand2",
            font=("Segoe UI", 10, "bold"),
            state="disabled"
        )
        self.send_btn.pack(side="left", padx=(6, 0))

        # ── строка статуса ────────────────────────────────────────────────
        self.status_var = tk.StringVar(value="Initialising…")
        status_bar = tk.Label(
            self, textvariable=self.status_var,
            bg="#181825", fg="#a6adc8",
            font=("Segoe UI", 9), anchor="w", padx=10, pady=3
        )
        status_bar.pack(fill="x", side="bottom")

    # ── IO redirect ─────────────────────────────────────────────────────────

    def _setup_io(self):
        self._gui_out = GuiOutput(self.console)
        self._gui_in  = GuiInput(
            self.console, self.entry, self.send_btn, self._input_queue
        )

    # ── Actions ──────────────────────────────────────────────────────────────

    def _send_input(self, _event=None):
        value = self.entry.get()
        self.entry.delete(0, tk.END)
        self._input_queue.put(value)

    def _browse_file(self):
        """Открыть файл и вставить путь в строку ввода."""
        path = filedialog.askopenfilename(
            title="Select PGS score file",
            filetypes=[
                ("All supported", "*.txt *.csv *.tsv *.gz"),
                ("Text files", "*.txt"),
                ("CSV files", "*.csv"),
                ("TSV files", "*.tsv"),
                ("Gzipped", "*.gz"),
                ("All files", "*.*"),
            ]
        )
        if path:
            # Если поле ввода активно — вставляем туда
            if self.entry["state"] == "normal":
                self.entry.delete(0, tk.END)
                self.entry.insert(0, path)
            else:
                # иначе копируем в буфер обмена
                self.clipboard_clear()
                self.clipboard_append(path)
                messagebox.showinfo(
                    "Path copied",
                    f"Path copied to clipboard:\n{path}\n\n"
                    "Paste it when the script asks for the file name."
                )

    def _run_script(self):
        """Запускает main() из pgs_to_plink в фоновом потоке."""
        self.status_var.set("⚙️  Running…")
        self._log("=" * 60 + "\n")
        self._log("  PGS → PLINK Converter GUI  |  pgs_to_plink.py v5.0\n")
        self._log("=" * 60 + "\n\n")

        def target():
            old_stdout, old_stderr, old_stdin = sys.stdout, sys.stderr, sys.stdin
            sys.stdout = self._gui_out
            sys.stderr = self._gui_out
            sys.stdin  = self._gui_in
            try:
                # импортируем и запускаем main() из основного скрипта
                import importlib.util
                script_path = os.path.join(SCRIPT_DIR, "pgs_to_plink.py")
                spec = importlib.util.spec_from_file_location("pgs_to_plink", script_path)
                mod  = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(mod)
                # exec_module уже запустит if __name__ == "__main__" блок?
                # Нет — запускаем main() явно
                mod.main()
            except SystemExit:
                pass
            except Exception as exc:
                sys.stdout.write(f"\n❌ Unexpected error: {exc}\n")
            finally:
                sys.stdout, sys.stderr, sys.stdin = old_stdout, old_stderr, old_stdin
                self.after(0, self._on_done)

        self._script_thread = threading.Thread(target=target, daemon=True)
        self._script_thread.start()

    def _on_done(self):
        self.status_var.set("✅  Done. Press 'Restart' to run again.")
        # разблокировать поле ввода на случай, если скрипт завершился неожиданно
        self.entry.configure(state="disabled")
        self.send_btn.configure(state="disabled")

    def _restart(self):
        if self._script_thread and self._script_thread.is_alive():
            if not messagebox.askyesno(
                "Restart", "Script is still running. Force restart?"
            ):
                return
        # очищаем очередь
        while not self._input_queue.empty():
            try:
                self._input_queue.get_nowait()
            except queue.Empty:
                break
        self.console.configure(state="normal")
        self.console.delete("1.0", tk.END)
        self.console.configure(state="disabled")
        self._run_script()

    # ── Helpers ──────────────────────────────────────────────────────────────

    def _log(self, text: str):
        self.console.configure(state="normal")
        self.console.insert(tk.END, text)
        self.console.see(tk.END)
        self.console.configure(state="disabled")


# ─── Entry point ─────────────────────────────────────────────────────────────

if __name__ == "__main__":
    app = PgsGuiApp()
    app.mainloop()
